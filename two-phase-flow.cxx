

#include <bout/physicsmodel.hxx>
#include <derivs.hxx>
#include <invert_laplace.hxx>
#include <field_factory.hxx>
#include <bout/constants.hxx>

const Field3D Div_VOF(const Field3D &n, const Field3D &f, BoutReal D=1.0);
const Field3D Curvature(const Field3D &c, BoutReal eps=1e-5);
const Field3D Curvature_HF(const Field3D &c, BoutReal eps=1e-5);
const Field3D smoothXZ(const Field3D &f);
const Field3D Div_Perp_Lap_FV(const Field3D &a, const Field3D &f, bool xflux);

/// Simulates incompressible flow of two fluids
/// 
class TwoPhase : public PhysicsModel {
private:
  Field3D vorticity; ///< Curl(momentum) with v = Curl(psi)
  Field3D vof, vof2;   ///< Volume of Fluid for fluid 1 and fluid 2.
  // The fraction of fluid 0 is (1-vof-vof2)
  
  Field3D psi; ///< Stream function, calculated from vorticity
  
  Field3D density; ///< Mass density. Calculated from vof at each iteration
  Field3D viscosity; ///< Kinematic viscosity. Calculated from vof
  
  Field3D kappa; ///< Surface curvature
  Field3D surf_force; ///< Surface forcing term in vorticity
  int curv_method; ///< Curvature calculation method
  int visc_method; ///< Viscosity method
  
  BoutReal density0, density1, density2; ///< Density of each fluid
  FieldGeneratorPtr viscosity0_generator, viscosity1_generator, viscosity2_generator;
  
  BoutReal surface_tension; ///< N/m. Note: Only surface tension between fluid 0 and fluid 1

  BoutReal gravity; ///< Acceleration due to gravity
  
  Laplacian *laplace; // Laplacian inversion to get stream function
  
  bool boussinesq; ///< Assume constant density in inertia?
  BoutReal vof_D; ///< Anti-diffusion for VOF advection

  Field3D Svort;  ///< External vorticity source
  FieldGeneratorPtr Svort_generator; ///< Calculates the source in time
  bool have_Svort; /// True if there is an external source

  Field3D shear_rate; ///< Shear rate (gamma) in s^-1
  
  // No-slip boundary conditions
  bool no_slip_x_down, no_slip_x_up;
protected:
  /// Initialise simulation
  ///
  /// @param[in] restarting  Is this simulation restarting?
  ///
  /// @returns zero on success
  int init(bool UNUSED(restarting)) {
    // Read input options
    auto opt = Options::root()["model"];
    OPTION(opt, density0, 1.0); // Water
    OPTION(opt, density1, 0.1); // Air
    OPTION(opt, density2, 10.0); // Sand

    // Generators so that the viscosity can be a function of time
    viscosity0_generator =
      FieldFactory::get()->parse(opt["viscosity0"].withDefault<std::string>("1.0"), &opt);
    viscosity1_generator =
        FieldFactory::get()->parse(opt["viscosity1"].withDefault<std::string>("0.1"), &opt);
    viscosity2_generator =
        FieldFactory::get()->parse(opt["viscosity2"].withDefault<std::string>("0.1"), &opt);
    
    OPTION(opt, surface_tension, 0.0);
    OPTION(opt, curv_method, 0);
    OPTION(opt, visc_method, 1);

    OPTION(opt, no_slip_x_down, false);
    OPTION(opt, no_slip_x_up, false);
    
    OPTION(opt, gravity, 0.1);
   
    OPTION(opt, boussinesq, true);
    OPTION(opt, vof_D, 0.1);
    
    // Specify evolving variables
    SOLVE_FOR(vorticity, vof, vof2);
    
    // Save the stream function at each output
    SAVE_REPEAT(psi);
    
    // Save the curvature and surface forcing term
    SAVE_REPEAT(kappa, surf_force);

    SAVE_REPEAT(viscosity);
    SAVE_REPEAT(shear_rate);
    
    // Create Laplacian inversion solver
    laplace = Laplacian::create(&Options::root()["laplace"]);
    
    // Allocate memory for density and viscosity
    // since we set by index in the rhs() function
    density.allocate();
    viscosity.allocate();

    // Make sure vof is between 0 and 1
    for (const auto &i : vof) {
      if (vof[i] < 0.0) {
        vof[i] = 0.0;
      }else if (vof[i] > 1.0) {
        vof[i] = 1.0;
      }
      if (vof2[i] < 0.0) {
        vof2[i] = 0.0;
      }else if (vof2[i] > 1.0) {
        vof2[i] = 1.0;
      }
    }

    // Check if there is a source of vorticity
    have_Svort = Options::root()["vorticity"]["source"].isSet();
    if (have_Svort) {
      // Read vorticity source from input
      Svort_generator =
          FieldFactory::get()->parse(Options::root()["vorticity"]["source"],
                                     &Options::root()["vorticity"]);
      Svort.allocate(); // Values will be set in rhs()

      SAVE_REPEAT(Svort);
    } else {
      Svort = 0.0;
    }

    psi = 0.0; // Starting value for Picard
    
    return 0;
  }

  /// Calculate time derivatives of evolving variables
  ///
  /// @param[in] time  The simulation time
  ///
  /// @returns zero on success
  int rhs(BoutReal time) {
    // Communicate guard cells
    mesh->communicate(vof, vof2, vorticity);

    Coordinates *coord = mesh->getCoordinates();
    
    if (have_Svort) {
      // Calculate vorticity source at this time
      BOUT_FOR(i, Svort.getRegion("RGN_ALL")) {
        Svort[i] =
            Svort_generator->generate(mesh->GlobalX(i.x()), TWOPI * mesh->GlobalY(i.y()),
                          TWOPI * static_cast<BoutReal>(i.z()) /
                              static_cast<BoutReal>(mesh->LocalNz),
                          time);
      }
    }
    
    // Calculate density and viscosity, given vof
    // Need to calculate fluid fractions
    // 
    // Note: Dependency density -> psi -> shear_rate -> viscosity
    Field3D fluid0, fluid1, fluid2; 
    fluid0.allocate();
    fluid1.allocate();
    fluid2.allocate();
    for (const auto &i : vof) {
      fluid1[i] = vof[i];
      // Make sure fraction of fluid is between 0 and 1
      if (fluid1[i] < 0.0) {
        fluid1[i] = 0.0;
      } else if (fluid1[i] > 1.0) {
        fluid1[i] = 1.0;
      }

      fluid2[i] = vof2[i];
      // Make sure fraction of fluid is between 0 and 1
      if (fluid2[i] < 0.0) {
        fluid2[i] = 0.0;
      } else if (fluid2[i] > 1.0) {
        fluid2[i] = 1.0;
      }

      // Check that the sum is <= 1
      BoutReal sum12 = fluid1[i] + fluid2[i];
      if (sum12 > 1.0) {
        fluid1[i] /= sum12;
        fluid2[i] /= sum12;
      }

      // Fluid0 is whatever fraction is left
      fluid0[i] = 1.0 - fluid1[i] - fluid2[i];

      // Density and viscosity now a weighted sum of contributions from each fluid
      density[i] = fluid0[i] * density0 + fluid1[i] * density1 + fluid2[i] * density2;
    }
    
    // Calculate curvature
    if (curv_method == 0) {
      // Finite Differences, using VOF function
      // Note we have to smooth the VOF function
      // to reduce spurious curvatures
    
      Field3D vof_smooth = vof;
      for (int i=0;i<6;i++) {
        vof_smooth = smoothXZ(vof_smooth);
      }
      kappa = Curvature(vof_smooth);
    } else if (curv_method == 1) {
      // Use Height Function method
      kappa = Curvature_HF(vof);
    } else {
      throw BoutException("Unrecognised curvature calculation method");
    }
    mesh->communicate(kappa);

    // Calculate stream function
    if (boussinesq) {
      // Ignore density variations in inertia
      psi = laplace->solve(vorticity / density0);
    } else {
      // Picard iteration
      // Vort = Div( density * Grad(psi))
      //      = Div(density0 * Grad(psi)) + Div( (density - density0) * Grad(psi))
      // so
      // Div(density0 * Grad(psi)) = Vort - Div( (density - density0) * Grad(psi))
      // Since density0 is a constant, the LHS can be solved using the
      // same solver as the Boussinesq approximation.
      // The RHS is calculated using the last value of psi, and iterated
      // to find a self-consistent solution
      
      // Record the maximum value of psi at this iteration
      // and the last iteration
      BoutReal psi_max_old, psi_max = max(psi, true);
      
      BoutReal absdiff; // The absolute difference between psi_max and psi_max_old
      do {
        mesh->communicate(psi);
        Field3D rhs = vorticity - Div_Perp_Lap_FV((density - density0), psi, false);
        psi = laplace->solve(rhs / density0);

        psi_max_old = psi_max;
        psi_max = max(psi, true);
        absdiff = abs(psi_max_old - psi_max);
      } while( (absdiff/psi_max > 1e-5) && (absdiff > 1e-9));
    }
    mesh->communicate(psi); // Communicates guard cells

    // Calculate shear rate (gamma)
    // 
    // shear_rate = sqrt( (psi_xx - psi_zz)^2 + psi_xz^2 )
    //
    // Based on http://dx.doi.org/10.1615/AtoZ.n.non-newtonian_fluids
    shear_rate.allocate();
    BOUT_FOR(i, shear_rate.getRegion("RGN_NOBNDRY")) {
      // Mesh spacing
      BoutReal dx = coord->dx[i];
      BoutReal dz = coord->dz;

      // Offset indices
      auto ixp = i.xp(); // X+1
      auto ixm = i.xm(); // X-1
      auto izp = i.zp(); // Z+1
      auto izm = i.zm(); // Z-1

      // Second derivative in X
      BoutReal psi_xx = (psi[i.xp()] - 2. * psi[i] + psi[i.xm()]) / SQ(dx);
      // Second derivative in Z
      BoutReal psi_zz = (psi[izp] - 2. * psi[i] + psi[izm]) / SQ(dz);
      // Mixed derivative in X-Z
      BoutReal psi_xz =
          ((psi[ixp.zp()] - psi[ixm.zp()]) - (psi[ixp.zm()] - psi[ixm.zm()])) /
          (4. * dx * dz);

      shear_rate[i] = sqrt(2.*(SQ(psi_xx - psi_zz) + SQ(psi_xz)));
      
      // Now calculate viscosity, which may use shear rate

      // Update kinematic viscosity coefficients at the current time
      // Note: The variable "y" is (ab)used to represent the shear rate
      BoutReal viscosity0 = viscosity0_generator->generate(0, shear_rate[i], 0, time);
      BoutReal viscosity1 = viscosity1_generator->generate(0, shear_rate[i], 0, time);
      BoutReal viscosity2 = viscosity2_generator->generate(0, shear_rate[i], 0, time);
      
      /// Kinematic viscosity
      viscosity[i] = fluid0[i] * viscosity0 + fluid1[i] * viscosity1 + fluid2[i] * viscosity2;
    }
    
    // Boundary condition on vorticity
    // Thom's formula, imposing no-slip conditions on X boundaries
    // https://web.math.princeton.edu/~weinan/papers/cfd5.pdf
    if (no_slip_x_down && mesh->firstX()) {
      // Loop over the boundary at x = xstart -1/2
      for (int y = mesh->ystart; y <= mesh->yend; y++) {
        for (int z = 0; z < mesh->LocalNz; z++) {
          // Value of vorticity on the boundary
          BoutReal vort_bndry = density(mesh->xstart, y, z) * 2 *
                                psi(mesh->xstart, y, z) /
                                SQ(0.5 * coord->dx(mesh->xstart, y));

          vorticity(mesh->xstart - 1, y, z) =
              2 * vort_bndry - vorticity(mesh->xstart, y, z);
        }
      }
    }

    if (no_slip_x_up && mesh->lastX()) {
      // Loop over the boundary at x = xend +1/2
      for (int y = mesh->ystart; y <= mesh->yend; y++) {
        for (int z = 0; z < mesh->LocalNz; z++) {
          // Value of vorticity on the boundary
          // Note that spacing between boundary and cells is dx/2
          BoutReal vort_bndry = density(mesh->xend, y, z) * 2 *
                                psi(mesh->xend, y, z) /
                                SQ(0.5 * coord->dx(mesh->xend, y));

          vorticity(mesh->xend + 1, y, z) =
              2 * vort_bndry - vorticity(mesh->xend, y, z);
        }
      }
    }

    
    // Vof, advected by flow
    ddt(vof) = -Div_VOF(vof, psi, vof_D);
    ddt(vof2) = -Div_VOF(vof2, psi, vof_D);
    
    // vorticity equation
    surf_force = surface_tension * bracket(kappa, vof, BRACKET_ARAKAWA);

    ddt(vorticity) =
      - bracket(psi, vorticity, BRACKET_ARAKAWA)
      - gravity * DDZ(density)
      - surf_force
      + Svort
      ;

    
    if (visc_method == 0) {
      // Simple viscosity term, assuming constant viscosity
      ddt(vorticity) += viscosity*Delp2(vorticity);
      
    } else if (visc_method == 1) {
      // Corrected form of viscosity. Adds a term which should cancel
      // the energy-violation due to shearing boundaries
      
      ddt(vorticity) +=  Div_Perp_Lap_FV(viscosity, vorticity, false);
      
    } else {
      // Full form of viscosity tensor, including variations in mu
      // (density or viscosity).
      // NOTE: Doesn't seem to be well behaved numerically.

      Field3D mu = viscosity * density;
      
      Field3D delp2_psi = Delp2(psi);
      Field3D ddz_psi = DDZ(psi);
      Field3D ddx_psi = DDX(psi);
      Field3D ddz_mu = DDZ(mu);
      Field3D ddx_mu = DDX(mu);

      mesh->communicate(delp2_psi, ddz_psi, ddx_psi, ddz_mu, ddx_mu);
      delp2_psi.applyBoundary("neumann");
      ddz_psi.applyBoundary("neumann");
      ddx_psi.applyBoundary("neumann");
      ddz_mu.applyBoundary("neumann");
      ddx_mu.applyBoundary("neumann");

      ddt(vorticity) += Delp2(mu * delp2_psi)
        + bracket(ddx_mu, ddz_psi, BRACKET_ARAKAWA)
        + bracket(ddx_psi, ddz_mu, BRACKET_ARAKAWA)
        ;
    }
    
    return 0;
  }
};


BoutReal minmod_vof(BoutReal a, BoutReal b, BoutReal c) {
  BoutReal sb = SIGN(b);
  return sb * BOUTMAX(0.0, BOUTMIN(sb*a, fabs(b), sb*c));
}

/// Smoothing in X-Z plane, using 9 point stencil
const Field3D smoothXZ(const Field3D &f) {
  Field3D result;
  
  result = copy(f);
  
  for(const auto &i : result.getRegion(RGN_NOBNDRY)) {
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

  Coordinates *coord = mesh->getCoordinates();
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
          
          // Need to account for:
          // 1) The velocity, since this flux is counteracting
          //    the diffusion caused by the advection
          // 2) The surface normal direction. The flux should
          //    be in the direction normal to the surface, and
          //    minimise distortion of the curvature of the surface
          //    since this will lead to parasitic flows

#if 1
          BoutReal coef = D*fabs(vL);

          // Calculate normal vector at this face
          BoutReal nx = (c - m)/dx(i,j);
          BoutReal nz = 0.25*(n(i,j,kp) + n(i-1,j,kp) - n(i,j,km) - n(i-1,j,km)) / dz;

          BoutReal nmag = sqrt(SQ(nx) + SQ(nz));
          if (nmag > 1e-6) {
            // Project flux onto normal vector
            coef *= fabs(nx)/nmag;
          }
#else
          BoutReal coef = D;
#endif

          flux += coef*minmod_vof(p - c, c - m, m - mm);

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

#if 1
        BoutReal coef = D*fabs(vD);

        BoutReal nz = (c - m)/dz;
        BoutReal nx = 0.25*(n(i+1,j,k) + n(i+1,j,km) - n(i-1,j,k) - n(i-1,j,km))/dx(i,j);
        BoutReal nmag = sqrt(SQ(nx) + SQ(nz));
        if (nmag > 1e-6) {
          // Project flux onto normal vector
          coef *= fabs(nz)/nmag;
        }
#else
        BoutReal coef = D;
#endif

        flux += coef*minmod_vof(p - c, c - m, m - mm);
        
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
  
  Coordinates *coord = mesh->getCoordinates();
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

        // Can't resolve radius of curvature smaller than
        // the grid size, so limit the magnitude of the
        // curvature
        BoutReal kmax = 0.5*(1./dz + 1./dx(i,j));
        if (result(i,j,k) < -kmax) {
          result(i,j,k) = -kmax;
        }
        if(result(i,j,k) > kmax) {
          result(i,j,k) = kmax;
        }

      }
  return result;
}

// Height Function estimation of the curvature
// From M.Sussman JCP 187 (2003) 110-136
// "A second order coupled level set and volume-of-fluid
// method for computing growth and collapse of vapor bubbles"
//
const Field3D Curvature_HF(const Field3D &c, BoutReal eps) {

  Field3D result(0.0);

  Coordinates *coord = mesh->getCoordinates();
  Field2D J = coord->J;
  Field2D dx = coord->dx;
  BoutReal dz = coord->dz;
  
  for(int i=mesh->xstart;i<=mesh->xend+1;i++)
    for(int j=mesh->ystart;j<=mesh->yend;j++)
      for(int k=0;k<mesh->LocalNz;k++) {
        
        int kp = (k+1) % (mesh->LocalNz);
        int km = (k-1+mesh->LocalNz) % (mesh->LocalNz);
        
        // Determine orientation of boundary
        
        /*
        BoutReal nx = (c(i+1,j,km) + c(i+1,j,k) + c(i+1,j,kp) 
                       - c(i-1,j,km) - c(i-1,j,k) - c(i-1,j,kp));
        BoutReal nz = ( c(i-1,j,kp) + c(i,j,kp) + c(i+1,j,kp)
                        - c(i-1,j,km) - c(i,j,km) - c(i+1,j,km) );
        */
        
        BoutReal nx = (c(i+1,j,k) - c(i-1,j,k))/(2.*dx(i,j,k));
        BoutReal nz = (c(i,j,kp) - c(i,j,km)) / (2.*dz);
        
        BoutReal hm, hc, hp; // Heights at -(m), centre (c) and +(p)
        BoutReal d;
        if(abs(nx) > abs(nz)) {
          // Difference in X bigger than difference in Z
          // so normal is mainly in X direction
          
          // Sum the height function
          hm = c(i,j,km);
          hc = c(i,j,k);
          hp = c(i,j,kp);
          for(int m=1;m<=3;m++) {
            if(i-m >= 0) {
              hm += c(i-m,j,km);
              hc += c(i-m,j,k);
              hp += c(i-m,j,kp);
            }
            if(i+m < mesh->LocalNz) {
              hm += c(i+m,j,km);
              hc += c(i+m,j,k);
              hp += c(i+m,j,kp);
            }
          }
          /*
          // This seemed like a good idea, but gives wrong result
          // for capillary wave test
          if (hc < 1.0) {
            continue;
          }
          */
          
          d = dz;
          
          hm *= dx(i,j);
          hc *= dx(i,j);
          hp *= dx(i,j);
          
        }else {
          // Normal mainly in Z direction
          
          hm = c(i-1,j,km) + c(i-1,j,k) + c(i-1,j,kp);
          hc = c(i,j,km) + c(i,j,k) + c(i,j,kp);
          hp = c(i+1,j,km) + c(i+1,j,k) + c(i+1,j,kp);
          
          for(int m=2;m<=3;m++) {
            kp = (kp+1) % mesh->LocalNz;
            km = (km-1+mesh->LocalNz) % mesh->LocalNz;
            
            hm += c(i-1,j,km) + c(i-1,j,kp);
            hc += c(i,j,km) + c(i,j,kp);
            hp += c(i+1,j,km) + c(i+1,j,kp);
          }
          
          /*
          if (hc < 1.0) {
            continue;
          }
          */

          hm *= dz;
          hc *= dz;
          hp *= dz;
          
          d = dx(i,j);
        }

        // First derivative of height
        BoutReal dh = (hp - hm) / (2.*d);
        
        // Second derivative of height
        BoutReal d2h = (hp - 2.*hc + hm)/SQ(d); 
        
        // Curvature of line in 2D
        result(i,j,k) = d2h / pow(1.0 + SQ(dh), 3./2);

        // Can't resolve radius of curvature smaller than
        // the grid size, so limit the magnitude of the
        // curvature
        BoutReal kmax = 1./d;
        if (result(i,j,k) < -kmax) {
          result(i,j,k) = -kmax;
        }
        if(result(i,j,k) > kmax) {
          result(i,j,k) = kmax;
        }
      }
  return result;
}

const Field3D Div_Perp_Lap_FV(const Field3D &a, const Field3D &f, bool xflux) {

  Field3D result = 0.0;

  Coordinates *coord = mesh->getCoordinates();
  
  //////////////////////////////////////////
  // X-Z diffusion
  //
  //            Z
  //            |
  //
  //     o --- gU --- o
  //     |     nU     |
  //     |            |
  //    gL nL      nR gR    -> X
  //     |            |
  //     |     nD     |
  //     o --- gD --- o
  //

  Field3D fs = f;
  Field3D as = a;

  for (int i = mesh->xstart; i <= mesh->xend; i++)
    for (int j = mesh->ystart; j <= mesh->yend; j++)
      for (int k = 0; k < mesh->LocalNz; k++) {
        int kp = (k + 1) % mesh->LocalNz;
        int km = (k - 1 + mesh->LocalNz) % mesh->LocalNz;

        // Calculate gradients on cell faces

        BoutReal gR = (coord->g11(i, j) + coord->g11(i + 1, j)) *
                          (fs(i + 1, j, k) - fs(i, j, k)) /
                          (coord->dx(i + 1, j) + coord->dx(i, j)) +
                      0.5 * (coord->g13(i, j) + coord->g13(i + 1, j)) *
                          (fs(i + 1, j, kp) - fs(i + 1, j, km) + fs(i, j, kp) -
                           fs(i, j, km)) /
                          (4. * coord->dz);

        BoutReal gL = (coord->g11(i - 1, j) + coord->g11(i, j)) *
                          (fs(i, j, k) - fs(i - 1, j, k)) /
                          (coord->dx(i - 1, j) + coord->dx(i, j)) +
                      0.5 * (coord->g13(i - 1, j) + coord->g13(i, j)) *
                          (fs(i - 1, j, kp) - fs(i - 1, j, km) + f(i, j, kp) -
                           f(i, j, km)) /
                          (4. * coord->dz);

        BoutReal gD = coord->g13(i, j) * (fs(i + 1, j, km) - fs(i - 1, j, km) +
                                         fs(i + 1, j, k) - fs(i - 1, j, k)) /
                          (4. * coord->dx(i, j)) +
                      coord->g33(i, j) * (fs(i, j, k) - fs(i, j, km)) / coord->dz;

        BoutReal gU = coord->g13(i, j) * (fs(i + 1, j, kp) - fs(i - 1, j, kp) +
                                         fs(i + 1, j, k) - fs(i - 1, j, k)) /
                          (4. * coord->dx(i, j)) +
                      coord->g33(i, j) * (fs(i, j, kp) - fs(i, j, k)) / coord->dz;

        // Flow right
        BoutReal flux = gR * 0.25 * (coord->J(i + 1, j) + coord->J(i, j)) *
                        (as(i + 1, j, k) + as(i, j, k));

        result(i, j, k) += flux / (coord->dx(i, j) * coord->J(i, j));
        // result(i+1,j,k) -= flux / (coord->dx(i+1,j)*coord->J(i+1,j));

        // Flow left
        flux = gL * 0.25 * (coord->J(i - 1, j) + coord->J(i, j)) *
               (as(i - 1, j, k) + as(i, j, k));

        result(i, j, k) -= flux / (coord->dx(i, j) * coord->J(i, j));
        // result(i-1,j,k) += flux / (coord->dx(i+1,j)*coord->J(i+1,j));

        // Flow up

        flux = gU * 0.5 * (as(i, j, k) + as(i, j, kp)) / coord->dz;
        result(i, j, k) += flux;
        // result(i,j,kp) -= flux;

        flux = gD * 0.5 * (as(i, j, k) + as(i, j, km)) / coord->dz;
        result(i, j, k) -= flux;
        // result(i,j,km) += flux;
      }
 
  return result;
}

BOUTMAIN(TwoPhase);
