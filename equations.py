from scipy.optimize import fsolve
import mpmath


def solve(gamma,k_bar_obs,k_max_obs) :
    def equations(p):
        alpha, kappa0,kappaC,P0,k_bar = p
        return (alpha-1+(kappa0/kappaC)**(gamma-2),
            P0- (gamma-1)*(alpha*kappa0)**(gamma-1)*mpmath.gammainc(1-gamma,alpha*kappa0),
            k_max_obs-alpha*kappaC,
            k_bar*alpha**2-(1-P0)*k_bar_obs,
            kappa0*(gamma-1)-k_bar*(gamma-2))
    return  fsolve(equations, (1, 1,k_max_obs,0.5,k_bar_obs))

def verify(alpha, kappa0,kappaC,P0,k_bar,gamma,k_bar_obs,k_max_obs) :
    def equations(p):
        alpha, kappa0,kappaC,P0,k_bar = p
        return (alpha-1+(kappa0/kappaC)**(gamma-2),
            P0- (gamma-1)*(alpha*kappa0)**(gamma-1)*mpmath.gammainc(1-gamma,alpha*kappa0),
            k_max_obs-alpha*kappaC,
            k_bar*alpha**2-(1-P0)*k_bar_obs,
            kappa0*(gamma-1)-k_bar*(gamma-2))
    print equations((alpha, kappa0,kappaC,P0,k_bar))
