import powerlaw
import numpy as np
import matplotlib.pyplot as plt

def get_exponent(graphe):
    return 1+(graphe.number_of_nodes()/sum(np.log(graphe.degree().values())))

def plot_gamma(graphe) :
    #fitness = powerlaw.Fit(graphe.degree().values(),xmin =1)
    #print fitness.alpha,fitness.xmin
    alpha = 1+(graphe.number_of_nodes()/sum(np.log(graphe.degree().values())))

    #powerlaw.plot_pdf(graphe.degree().values())
    xdata,ydata = powerlaw.pdf(graphe.degree().values())

    vala ={}
    xvala=[]
    yvala=[]
    somme = 0
    exceptions =[]
    last_exception = -1
    for value in graphe.degree().values() :
        somme+=1
        try :
            vala[value]+=1
        except :
            vala[value] =1

    for degree in range(1,max(vala.keys())+1) :
        try :
            yvala.append(float(vala[degree])/somme)
            xvala.append(degree)

        except :
            if len(exceptions) ==0 :
                last_exception =degree
                #print degree
            exceptions.append(degree)
    #slope, intercept, r_value, p_value, std_err = stats.linregress(
        #np.log(xvala[:last_exception]),np.log(yvala[:last_exception]))

    #print slope

    x = xvala

    plt.plot(x,np.power(x,-alpha),color = "orange")
    #plt.plot(x,np.power(x,-2.34296138661),color ="blue")
    #plfit.plfit(graphe.degree().values(),quiet = False,verbose = True,silent = False)
    plt.plot(xvala,yvala,color = "black")
    plt.yscale('log')
    plt.xscale('log')
    #plt.ylim([0,0.4])
    #plt.xlim([0,20])
    plt.show()
    #xdata,ydata=powerlaw.pdf(graphe.degree().values())