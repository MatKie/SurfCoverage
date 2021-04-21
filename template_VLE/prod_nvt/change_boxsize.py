import subprocess
from mkutils import PlotGromacs
import os

stats = PlotGromacs.get_gmx_stats(os.path.join('..', 'eq_npt', 'energies.out'))
box_z = stats.get('Box-Z').get('Average')
box_x = stats.get('Box-X').get('Average')

process = subprocess.run(['gmx editconf -f ../eq_npt/confout.gro -o input.gro \
                          -box {0:.5f} {0:.5f} {1:.5f}'.format(box_x, box_z)],
                        check=True, 
                        shell=True
                        )

