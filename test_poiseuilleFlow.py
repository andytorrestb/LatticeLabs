import BoundaryConditions as BC
import Simulation as sim
import D2Q9
import ComputationalDomain as CD
import PostProcessor

if __name__ == "__main__":
    # Initialize the D2Q9 lattice
    tau = 1.0
    d2q9 = D2Q9.D2Q9(tau)

    # Initialize the computational domain as a Poiseuille flow.
    nx = 100 
    ny = int(nx / 2)
    dx = 1.0
    dy = 1.0
    cd = CD.PoiseuilleFlow(nx, ny, dx, dy)
    cd.set_lattice(d2q9)
    cd.set_inlet_density(5.0)
    cd.initialize()

    # Initialize the simulation
    simulation = sim.Simulation(d2q9, cd)
    results = simulation.run(100)

    print(results)

    # Post-process the results
    pp = PostProcessor.PostProcessor(results, simulation)
    pp.process()

