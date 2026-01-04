import numpy as np
import matplotlib.pyplot as plt

def Simula():
    dt = 0.001 
    Periodo = 100.0
    VecTempi = np.arange(0, Periodo, dt)
    Passi = len(VecTempi)
    VecTheta = np.zeros(Passi)
    VecOmega = np.zeros(Passi)
    VecAcc = np.zeros(Passi)
    
    mu0 = 1.256e-6      
    massa = 0.0085          
    lung = 0.08           
    theta = 0.1        
    thetaPrima = 0.0   
    delta = 0.5        
    omega0 = 2.0       
    m = 0.0212        
    q = 15.0          
    a = 0.32     
    g = 9.8067 
    
    VecTheta[0] = theta
    VecOmega[0] = thetaPrima
    
    Delay = int(delta / dt)

    for i in range(Passi - 1):
        if i >= Delay:
            SenMemoria = np.sin(VecTheta[i - Delay])
        else:
            SenMemoria = 0.0 

        def Accelerazione(th, om):
            n1 = (3 * mu0) * (m * (np.sin(th) + SenMemoria))
            d1 = 2 * np.pi * (massa * np.sin(th) * lung * np.cos(th))
            p1 = n1 / d1
            d2 = (a + (lung * SenMemoria) + np.sin(th))**4
            p2 = 1 / d2
            
            n3 = (omega0 * massa) * (om * lung) * np.tan(th)
            d3 = q * massa
            p3 = n3 / d3

            gravita = -(g / lung) * np.sin(th)
            return (p1 * p2) - p3 + gravita

        # RK - 4 
        ThAtt = VecTheta[i]
        OmAtt = VecOmega[i]

        k1_v = Accelerazione(ThAtt, OmAtt)
        k1_x = OmAtt
        
        k2_v = Accelerazione(ThAtt + 0.5 * dt * k1_x, OmAtt + 0.5 * dt * k1_v)
        k2_x = OmAtt + 0.5 * dt * k1_v
        
        k3_v = Accelerazione(ThAtt + 0.5 * dt * k2_x, OmAtt + 0.5 * dt * k2_v)
        k3_x = OmAtt + 0.5 * dt * k2_v
        
        k4_v = Accelerazione(ThAtt + dt * k3_x, OmAtt + dt * k3_v)
        k4_x = OmAtt + dt * k3_v
        
        VecOmega[i+1] = OmAtt + (dt/6.0)*(k1_v + 2 * k2_v + 2 * k3_v + k4_v)
        VecTheta[i+1] = ThAtt + (dt/6.0)*(k1_x + 2 * k2_x + 2 * k3_x + k4_x)
        VecAcc[i] = k1_v


    plt.figure(figsize=(12, 6))
    plt.plot(VecTempi, VecOmega, label=r"$\theta'$ (Velocità)", color='blue', linewidth=1)
    plt.plot(VecTempi, VecAcc, label=r"$\theta''$ (Accelerazione)", color='red', alpha=0.7, linewidth=1)
    plt.axhline(0, color='black', linewidth=1)
    plt.title("Analisi Pendolo Magnetico")
    plt.xlabel("Tempo [s]")
    plt.ylabel("Ampiezza [rad]")
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.legend()
    fig = plt.gcf()
    fig.canvas.manager.set_window_title('Analisi Grafica')
    plt.show()

if __name__ == "__main__":
    Simula()