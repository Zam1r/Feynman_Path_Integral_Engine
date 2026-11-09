from plot_convergence import main as convergence
from plot_diagnostics import main as diagnostics
from plot_landscape import main as landscape
from plot_milestone import main as milestone
from plot_paths import main as paths


if __name__ == "__main__":
    milestone()
    convergence()
    paths()
    diagnostics()
    landscape()
