import wavesolve
from scipy.optimize import fmin_cobyla

def objective(args):
    wavesolve.solve(args, 10)

def constr1(args):
    return args[0] + args[2]

def constr2(args):
    return args[0] + args[4]

def constr3(args):
    return args[0] + args[4]


def main():
    #objective([.1490,1.199,.899,1.18,-.077,.2250])
    fmin_cobyla(objective, [.100,1.0,.899,1.0,-.077,.2250], [constr1, constr2], rhoend=1e-7)



if __name__ == '__main__':
    main()
