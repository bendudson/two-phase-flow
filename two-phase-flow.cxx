

#include <bout/physicsmodel.hxx>
#include <derivs.hxx>
#include <invert_laplace.hxx>

/// Simulates incompressible flow of two fluids
/// 
class TwoPhase : public PhysicsModel {
private:
  Field3D vorticity;
  Field3D vof;   ///< Volume of Fluid

  Field3D psi; ///< Stream function, calculated from vorticity
  
  Field3D density; ///< Mass density. Calculated from vof at each iteration
  Field3D viscosity; ///< Kinematic viscosity. Calculated from vof
  
  BoutReal density0, density1; // Density of each fluid
  BoutReal viscosity0, viscosity1; // Kinematic viscosity of each fluid

  BoutReal gravity;
  
  Laplacian *laplace; // Laplacian inversion to get stream function
  
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

    OPTION(opt, gravity, 0.1);
    
    // Specity evolving variables
    SOLVE_FOR2(vorticity, vof);

    // Save the stream function at each output
    SAVE_REPEAT(psi);
    
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
    
    // Calculate stream function
    psi = laplace->solve(vorticity / density);
    mesh->communicate(psi);
    
    // Vof, advected by flow
    ddt(vof) = -bracket(psi, vof, BRACKET_SIMPLE);

    // vorticity equation
    ddt(vorticity) =
      - bracket(psi, vorticity, BRACKET_ARAKAWA)
      - gravity * DDZ(density)
      + viscosity*Delp2(vorticity)
      ;
    return 0;
  }
};

BOUTMAIN(TwoPhase);
