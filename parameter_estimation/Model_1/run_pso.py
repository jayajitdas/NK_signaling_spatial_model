import os, sys, shutil, functools, hashlib, glob
import numpy as np
from time import sleep
from pyswarm import pso
from setup_full_run import setup_run
from std_density import calc_std_density
from S4_std_density import calc_S4_std_density
from pair_corr import calc_pair_corr
from S4_adjusted_pair_corr import S4_calc_pair_corr
from thres_S4_2_avg_voxel import S4_adjust_voxel_size
from voxel_size import adjust_voxel_size

debug_flag = False

def counted(f):
    @functools.wraps(f)
    def wrapped(*args, **kwargs):
        wrapped.calls += 1
        return f(*args, **kwargs)
    wrapped.calls = 0
    return wrapped

@counted
def myfunc(x):
    dir_name = "md5_"+hashlib.md5(str(x)).hexdigest()
    if debug_flag:
        print(dir_name)
        print(x)
        print('Attempting to load solution')
    # try loading in the answer from a previous run
    try:
        f = open(dir_name+"/sites.0.1",'r')
    except IOError:
        pass
    else:
        if debug_flag:
            print('Found solution')
        cost=0.0
        c_old=0.0
        sleep(5)
        f.close()

        sleep(10)
        f = open(dir_name+"/sites.0.1",'r')
        adjust_voxel_size(f,dir_name,0.0,60)
        f.close()

        name = "name"

        ff = open("sites.org_S4_1_60.60",'r')
        S4_adjust_voxel_size(ff,name,0.0,1)
        ff.close()

        f = open(dir_name+"/sites.0.0_"+str(dir_name)+"adjusted125.125_60",'r')
        height_mean,height_stdev = calc_std_density(f)
        f.close()

        f = open("sites.0.0_"+str(name)+"_thres_adjusted_S4_1_26.26",'r')
        height_mean_S4_2,height_stdev_S4_2 = calc_S4_std_density(f)
        f.close()

        cost = cost +  (( ( ((height_mean - height_mean_S4_2)/height_mean_S4_2) **2)/3.0) + ( (((height_stdev - height_stdev_S4_2)/height_stdev_S4_2) **2)/3.0) )

        popt1 = []
        f = open(dir_name+"/sites.0.0_"+str(dir_name)+"adjusted125.125_60",'r')
        popt1 = calc_pair_corr(f)
        f.close()
        
        image_popt1 = []
        ff = open("sites.0.0_"+str(name)+"_thres_adjusted_S4_1_26.26",'r')
        image_popt1 = S4_calc_pair_corr(ff)
        ff.close()
        diff_popt = 0.0
        avg_diff_popt = 0.0
        len_popt1 = min(len(popt1), len(image_popt1))
        print("sim",len(popt1),"image", len(image_popt1))
        sys.stdout.flush()
        for ii in range(len_popt1):
            diff_popt = (popt1[ii] - image_popt1[ii])**2
            avg_diff_popt = avg_diff_popt + diff_popt
        cost = cost + ( (avg_diff_popt)/3.0)
        
        if debug_flag:
            print(24,height_mean,height_CV,cost)
        
        print(dir_name,"function",height_mean,height_mean_S4_2,height_stdev,height_stdev_S4_2,avg_diff_popt,cost)
        sys.stdout.flush()
        print(dir_name,"x=%.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f"%tuple(x),"cost=%.8f"%cost)
        sys.stdout.flush()
        return cost

    if debug_flag:
        print('No old solution found. Attempting to start run')
    # if no previous answer exists, we may have to do a run
    try:
        # try making the directory for the new run
        os.mkdir(dir_name)
    except OSError:
        # if the directory already exists, then either a run is already in progress, or had an error

        # check if the run ended with an error by seeing if the slurm error file has anything in it
        for fname in glob.glob(dir_name+'/slurm*err'):
            if os.stat(fname).st_size > 0:
                if debug_flag:
                    print("Error in run. Returning 1e99")
                print(dir_name, "err")
                return 1e99
    else:
        # if the directory was created start the new run
        shutil.copy2("spk_redsky",dir_name)
        shutil.copy2("sites.30.30",dir_name)
        os.chdir(dir_name)
        setup_run(x)
        os.system("sbatch run.sh")
        sleep(2)
        os.chdir("..")
        if debug_flag:
            print('Started run')
        
    if debug_flag:
            print("returning None")
    print(dir_name, "running...")
    sys.stdout.flush()
    return None

np.random.seed(42)

lb = [-1, -1, -1, -2, -1, 0, -2, -1, -2, -2, -1, -1, -0.2, -0.1, -2, -0.69, -1]
ub = [+2, 0, +1, +1, +2, +2, +1, +2, +1, +1, +1, +1, +0.3, 0, +1, +0.5, +0.7]

xopt1, fopt1 = pso(myfunc, lb, ub,swarmsize=200,omega=0.5, phip=2.5, phig=1.5, maxiter=100, minstep=0.000001, minfunc=0.0000001, debug=True)

print('The optimum is at:')
print('    {}'.format(xopt1))
print('Optimal function value:')
print('    myfunc: {}'.format(fopt1))
print('    called {} times'.format(myfunc.calls))
