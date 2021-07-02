# Run a parameter sweep of the Held Suarez model
# by varying the rotation rate from 1% to 1000% of Earth's rot rate
import numpy as np
import os
from isca import Experiment, SocratesCodeBase, FailedRunError, GFDL_BASE
from isca.util import exp_progress

from socrates_aquaplanet import namelist, diag

cb = SocratesCodeBase.from_directory(GFDL_BASE)

namelist['main_nml'] = {
        'days'   : 30,
        'hours'  : 0,
        'minutes': 0,
        'seconds': 0,
        'dt_atmos':600,
        'current_date' : [1,1,1,0,0,0],
        'calendar' : 'thirty_day'
}

inputfiles = [os.path.join(GFDL_BASE,'input/rrtm_input_files/ozone_1990.nc'),
os.path.join(GFDL_BASE,'input/ideal_exo_0.nc'),
os.path.join(GFDL_BASE,'input/ideal_exo_1.nc'),
os.path.join(GFDL_BASE,'input/ideal_exo_2.nc'),
os.path.join(GFDL_BASE,'input/ideal_exo_3.nc'),
os.path.join(GFDL_BASE,'input/ideal_exo_4.nc'),
os.path.join(GFDL_BASE,'input/ideal_exo_5.nc'),
os.path.join(GFDL_BASE,'input/ideal_exo_6.nc'),
]



for s in range(7):
    exp_name = 'soc_topo_%d' % s
    exp = Experiment(exp_name, codebase=cb)
    exp.namelist = namelist.copy()
    exp.inputfiles = inputfiles
    exp.diag_table = diag
    if s > 0:
        exp.update_namelist({'spectral_init_cond_nml':{'topography_option':'input'}})
        exp.update_namelist({'spectral_init_cond_nml': {'topog_file_name': 'ideal_exo_%d.nc' % s}})
    try:
        # run with a progress bar with description showing omega
        with exp_progress(exp, description='o%.0f d{day}' % s) as pbar:
            exp.run(1, use_restart=False, num_cores=8)

        for n in range(2, 11):
            with exp_progress(exp, description='o%.0f d{day}' % s) as pbar:
                exp.run(n)
                exp.delete_restart(n-1)

    except FailedRunError as e:
        # don't let a crash get in the way of good science
        # (we could try and reduce timestep here if we wanted to be smarter)
        continue
