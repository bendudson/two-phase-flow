

#include <bout/physicsmodel.hxx>
#include <derivs.hxx>
#include <invert_laplace.hxx>

const Field3D Div_VOF(const Field3D &n, const Field3D &f, BoutReal D=1.0);
const Field3D Curvature(const Field3D &c, BoutReal eps=1e-5);
const Field3D smoothXZ(const Field3D &f);

/// Simulates incompressible flow of two fluids
/// 
class TwoPhase : public PhysicsModel {
private:
  Field3D vorticity; ///< Curl(momentum) with v = Curl(psi)
  Field3D vof;   ///< Volume of Fluid

  Field3D psi; ///< Stream function, calculated from vorticity
  
  Field3D density; ///< Mass density. Calculated from vof at each iteration
  Field3D viscosity; ///< Kinematic viscosity. Calculated from vof
  
  Field3D kappa; ///< Surface curvature
  Field3D surf_force; ///< Surface forcing term in vorticity
  Field3D vof_smooth; ///< Smoothed VOF used for curvature calculation

  BoutReal density0, density1; // Density of each fluid
  BoutReal viscosity0, viscosity1; // Kinematic viscosity of each fluid
  BoutReal surface_tension; // sigma

  BoutReal gravity;
  
  Laplacian *laplace; // Laplacian inversion to get stream function
  
  BoutReal vof_D; // Anti-diffusion for VOF advection
protected:
  /// Initialise simulation
  ///
  /// @param[in] restarting  Is this simulation restarting?
  ///
  /// @returns zero on success
  int init(bool restarting) {
    // Read input options
    Options *opt = Options::getRoot()->getSection("model");
    OPTION(opt, density0, 1.0);
    OPTION(opt, density1, 0.1);
    OPTION(opt, viscosity0, 1.0);
    OPTION(opt, viscosity1, 0.1);
    OPTION(opt, surface_tension, 0.0);

    OPTION(opt, gravity, 0.1);
    
    OPTION(opt, vof_D, 0.1);
    
    // Specity evolving variables
    SOLVE_FOR2(vorticity, vof);

    // Save the stream function at each output
    SAVE_REPEAT(psi);
    
    // Save the curvature
    SAVE_REPEAT3(kappa, surf_force, vof_smooth);

    // Create Laplacian inversion solver
    laplace = Laplacian::create(Options::getRoot()->getSection("laplace"));
    
    // Allocate memory for density and viscosity
    // since we set by index in the rhs() function
    density.allocate();
    viscosity.allocate();

    // Make sure vof is between 0 and 1
    for(const auto &i : vof) {
      if (vof[i] < 0.0) {
        vof[i] = 0.0;
      }else if (vof[i] > 1.0) {
        vof[i] = 1.0;
      }
    }
    return 0;
  }

  /// Calculate time derivatives of evolving variables
  ///
  /// @param[in] time  The simulation time
  ///
  /// @returns zero on success
  int rhs(BoutReal time) {
    // Communicate guard cells
    mesh->communicate(vof, vorticity);
    
    // Calculate density and viscosity, given vof
    
    for(const auto &i : vof) {
      BoutReal c = vof[i];
      // Make sure fraction of fluid is between 0 and 1
      if (c < 0.0) {
        c = 0.0;
      }else if (c > 1.0) {
        c = 1.0;
      }

      // When vof=0 -> density=density0 ; when vof=1 -> density=density1
      density[i] = c*density1 + (1.-c)*density0;
      viscosity[i] = c*viscosity1 + (1.-c)*viscosity0;
    }
    
    // Calculate curvature. Note we have to smooth the VOF function
    // to reduce spurious curvatures

    vof_smooth = vof;
    for(int i=0;i<6;i++) {
      vof_smooth = smoothXZ(vof_smooth);
    }
    kappa = Curvature(vof_smooth, 1e-9);
    mesh->communicate(kappa);

    // Calculate stream function
    psi = laplace->solve(vorticity / density);
    mesh->communicate(psi); // Communicates guard cells
    
    // Vof, advected by flow
    ddt(vof) = -Div_VOF(vof, psi, vof_D);

    // vorticity equation
    surf_force = surface_tension * bracket(kappa, vof, BRACKET_ARAKAWA);

    ddt(vorticity) =
      - bracket(psi, vorticity, BRACKET_ARAKAWA)
      - gravity * DDZ(density)
      + viscosity*Delp2(vorticity)
      - surf_force
      ;
    return 0;
  }
};


BoutReal minmod_vof(BoutReal a, BoutReal b, BoutReal c) {
  BoutReal sb = SIGN(b);
  return sb * BOUTMAX(0.0, BOUTMIN(sb*a, fabs(b), sb*c));
}

const Field3D smoothXZ(const Field3D &f) {
  Field3D result;
  
  result = copy(f);
  
  for(const auto &i : result.region(RGN_NOBNDRY)) {
    result[i] = 0.5*f[i] 
      + (
         f[i.xp()] + f[i.xm()]
         + f[i.zp()] + f[i.zm()] 
         + f[i.offset(1,0,1)] + f[i.offset(1,0,-1)]
         + f[i.offset(-1,0,1)] + f[i.offset(-1,0,-1)]
         ) / 16.;
  }
  mesh->communicate(result);
  return result;
}

/*
 * Divergence of a Volume Of Fluid (VOF) marker
 *  Div (n * b x Grad(f)/B)
 *
 * 
 * 
 */
const Field3D Div_VOF(const Field3D &n, const Field3D &f, BoutReal D) {
  Field3D result(0.0);
  
  //////////////////////////////////////////
  // X-Z advection.
  // 
  //             Z
  //             |
  // 
  //    fmp --- vU --- fpp
  //     |      nU      |
  //     |               |
  //    vL nL        nR vR    -> X
  //     |               |
  //     |      nD       |
  //    fmm --- vD --- fpm
  //

  Coordinates *coord = mesh->coordinates();
  Field2D J = coord->J;
  Field2D dx = coord->dx;
  BoutReal dz = coord->dz;
  
  for(int i=mesh->xstart;i<=mesh->xend+1;i++)
    for(int j=mesh->ystart;j<=mesh->yend;j++)
      for(int k=0;k<mesh->LocalNz;k++) {
        int kp = (k+1) % (mesh->LocalNz);
        int km = (k-1+mesh->LocalNz) % (mesh->LocalNz);
        int kmm = (km-1+mesh->LocalNz) % (mesh->LocalNz);
	
        // 1) Interpolate stream function f onto corners fmp, fpp, fpm
	
        BoutReal fmm = 0.25*(f(i,j,k) + f(i-1,j,k) + f(i,j,km) + f(i-1,j,km));
        BoutReal fmp = 0.25*(f(i,j,k) + f(i,j,kp) + f(i-1,j,k) + f(i-1,j,kp)); // 2nd order accurate
        BoutReal fpm = 0.25*(f(i,j,k) + f(i+1,j,k) + f(i,j,km) + f(i+1,j,km));
	
        // 2) Calculate velocities on cell faces vD and vL
        BoutReal vD = J(i,j)*(fmm - fpm)/dx(i,j); // -J*df/dx
        BoutReal vL = 0.5*(J(i,j)+J(i-1,j))*(fmp - fmm)/dz; // J*df/dz 

        // 3) Calculate n on the cell faces. The sign of the
        //    velocity determines which side is used.
	
        // X direction
	
        BoutReal c = n(i,j,k);
        BoutReal m = n(i-1,j,k);
        BoutReal mm = n(i-2,j,k);
        BoutReal p = n(i+1,j,k);
        
        // Left side
          
        if ((i==mesh->xstart) && (mesh->firstX())) {
          // At left boundary in X
          
        } else {
          // Not at a boundary
          
          // This implements 1st-order upwinding
          BoutReal flux;
          if (vL < 0.0) {
            // Flux from cell i into i-1
            flux = vL * c;
          } else {
            // Flux from cell i-1 into i
            flux = vL * m;
          }

          // Anti-diffusion. This algorithm from
          // "Anti-diffusion method for interface steepening in two-phase incompressible flow"
          // by K. K. So, X. Y. Hu and N. A. Adams
          
          flux += D*minmod_vof(p - c, c - m, m - mm);

          result(i,j,k)   -= flux / (dx(i,j) * J(i,j));
          result(i-1,j,k) += flux / (dx(i-1,j) * J(i-1,j));
        }
        
        m = n(i,j,km);
        mm = n(i,j,kmm);
        p = n(i,j,kp);
        
        BoutReal flux;
        if (vD < 0.0) {
          flux = vD * c;
        }else {
          flux = vD * m;
        }
        flux += D*minmod_vof(p - c, c - m, m - mm);
        
        flux /= J(i,j)*dz;
        
        result(i,j,k)   -= flux;
        result(i,j,km)  += flux;
      }
  return result;
}

/// Calculate the curvature
/// kappa = Div(n/|n|)
///
/// Calculation from: 
/// http://elie.korea.ac.kr/~cfdkim/papers/surface_tension.pdf
const Field3D Curvature(const Field3D &c, BoutReal eps) {
  Field3D result(0.0);
  
  Coordinates *coord = mesh->coordinates();
  Field2D J = coord->J;
  Field2D dx = coord->dx;
  BoutReal dz = coord->dz;
  
  for(int i=mesh->xstart;i<=mesh->xend+1;i++)
    for(int j=mesh->ystart;j<=mesh->yend;j++)
      for(int k=0;k<mesh->LocalNz;k++) {
        // Calculate normal vector at vertex of each cell
        // Need x and z components at (x,z) +/- (p/m) 1/2
        // 
        
        int kp = (k+1) % (mesh->LocalNz);
        int km = (k-1+mesh->LocalNz) % (mesh->LocalNz);

        BoutReal nx_pp = (c(i+1,j,k) + c(i+1,j,kp) - c(i,j,k) - c(i,j,kp))/(2.*dx(i,j));
        BoutReal nz_pp = (c(i,j,kp) + c(i+1,j,kp) - c(i,j,k) - c(i+1,j,k))/(2.*dz);
        
        BoutReal nx_mp = (c(i,j,k) + c(i,j,kp) - c(i-1,j,k) - c(i-1,j,kp))/(2.*dx(i-1,j));
        BoutReal nz_mp = (c(i-1,j,kp) + c(i,j,kp) - c(i-1,j,k) - c(i,j,k))/(2.*dz);

        BoutReal nx_pm = (c(i+1,j,km) + c(i+1,j,k) - c(i,j,km) - c(i,j,k))/(2.*dx(i,j));
        BoutReal nz_pm = (c(i,j,k) + c(i+1,j,k) - c(i,j,km) - c(i+1,j,km))/(2.*dz);

        BoutReal nx_mm = (c(i,j,km) + c(i,j,k) - c(i-1,j,km) - c(i-1,j,k))/(2.*dx(i,j));
        BoutReal nz_mm = (c(i-1,j,k) + c(i,j,k) - c(i-1,j,km) - c(i,j,km))/(2.*dz);
        
        // Get magnitudes

        BoutReal pp_mag = sqrt(SQ(nx_pp) + SQ(nz_pp));
        BoutReal pm_mag = sqrt(SQ(nx_pm) + SQ(nz_pm));
        BoutReal mp_mag = sqrt(SQ(nx_mp) + SQ(nz_mp));
        BoutReal mm_mag = sqrt(SQ(nx_mm) + SQ(nz_mm));

        // If any n vector is too small, set kappa to zero
        if ( (pp_mag < eps) || (pm_mag < eps) || (mp_mag < eps) || (mm_mag < eps)) {
            continue;
          }
        
        // Normalise
        
        nx_pp /= pp_mag;
        nz_pp /= pp_mag;
        
        nx_pm /= pm_mag;
        nz_pm /= pm_mag;
        
        nx_mp /= mp_mag;
        nz_mp /= mp_mag;

        nx_mm /= mm_mag;
        nz_mm /= mm_mag;

        // Calculate divergence

        result(i,j,k) = 0.5*(nx_pp + nx_pm - nx_mp - nx_mm)/dx(i,j)
          + 0.5*(nz_mp + nz_pp - nz_mm - nz_pm)/dz;

      }
  return result;
}

BOUTMAIN(TwoPhase);
