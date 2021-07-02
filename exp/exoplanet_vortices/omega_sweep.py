# Run a parameter sweep of the Held Suarez model
# by varying the rotation rate from 1% to 1000% of Earth's rot rate
import numpy as np
from isca import Experiment, DryCodeBase, FailedRunError, GFDL_BASE
from isca.util import exp_progress

from socrates_aquaplanet import namelist, diag

cb = DryCodeBase.from_directory(GFDL_BASE)

namelist['main_nml'] = {
        'days'   : 30,
        'hours'  : 0,
        'minutes': 0,
        'seconds': 0,
        'dt_atmos':600,
        'current_date' : [1,1,1,0,0,0],
        'calendar' : 'thirty_day'
}

earth_omega = 7.292e-5

scales = [1/16, 1/8, 1/4, 1/2, 1.0, 2.0, 4.0, 8.0, 16.0]

for s in scales:
    exp_name = 'soc_om_scale_%.0f' % s
    omega = earth_omega * (s/100.0)
    exp = Experiment(exp_name, codebase=cb)
    exp.namelist = namelist.copy()
    exp.diag_table = diag

    exp.update_namelist({'constants_nml': {'omega': omega}})
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
