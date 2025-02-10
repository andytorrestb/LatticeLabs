
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

class PostProcessor:
    def __init__(self, results, simulation):
        self.simulation = simulation


        self.results = results
        print(type(results))

        t = []
        u = []
        rho = []

        for result in results:
            t.append(result["time"])
            u.append(result["velocity"][0])
            rho.append(result["density"])

        self.t = np.array(t)
        self.u = np.array(u)
        self.rho = np.array(rho)
        return

    def create_meshgrid(self):
        x_min = 0
        x_max = self.simulation.computational_domain.nx
        nx = self.simulation.computational_domain.nx

        y_min = 0
        y_max = self.simulation.computational_domain.ny
        ny = self.simulation.computational_domain.ny
        
        x = np.linspace(x_min, x_max, nx)
        y = np.linspace(y_min, y_max, ny)

        X, Y = np.meshgrid(x, y)
        return X, Y

    def create_video(self, filename='output.mp4'):
        fig, ax = plt.subplots()
        line = ax.imshow(
            np.zeros(
                (self.simulation.computational_domain.ny, self.simulation.computational_domain.nx)
            ), 
            animated=True,
            label = 'Velocity',
            vmin = 0,
            vmax = 255
        )
        
        x_min = 0
        x_max = self.simulation.computational_domain.nx

        y_min = 0
        y_max = self.simulation.computational_domain.ny

        X, Y = self.create_meshgrid()

        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)
        ax.legend()

        def init():
            line.set_array(np.zeros((self.simulation.computational_domain.ny, self.simulation.computational_domain.nx)))
            return line,

        def update(frame):
            line.set_array(self.rho[frame])
            print(self.rho[frame])
            return line,

        ani = animation.FuncAnimation(fig, update, frames=len(self.t), init_func=init, blit=True, interval=30)
        ani.save(filename, writer='ffmpeg')
        
    def process(self):
        # do something with data
        self.create_video()
