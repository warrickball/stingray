import numpy as np
import scipy
from scipy import special
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from mpl_toolkits.mplot3d import axes3d
from scipy import signal
from scipy.fftpack import fft, fftshift
from scipy.fftpack import fft, fftshift
import sys
import accel_utils
from scipy.integrate import trapz
import time
#this function work just on chunks, never use it for FFT long more than 1e**3

def accelsearch( spectr, freq, T, vi = False, delta_rdot = 1 ):
    center = int( len( freq ) / 2 )
    # range rdot
    start_rdot = -20#-20*delta_rdot 
    end_rdot = 20#-start_rdot + delta_rdot
    range_rdot = np.arange(start_rdot,end_rdot, delta_rdot)
    print("min and max possible r_dot: " + str( delta_rdot / (T**2) )+ " " + str( np.max(range_rdot) / (T**2)))
    #maximum m index, 
    maximum_m_idx = 2 * int( np.max(range_rdot) )
    #centering the box in the interesting fequencys
    min_freq_idx = center - 2 * int ( maximum_m_idx / 2 )
    max_freq_idx = center + 2 * int ( maximum_m_idx / 2 ) + 1
    interest_freq = freq[ min_freq_idx : max_freq_idx ] 
    interest_spectr = spectr[ min_freq_idx : max_freq_idx ]
    #dimension where calculate the complex filter
    N = len( interest_freq )
    M = len( range_rdot )
    mesh = np.zeros( ( N, M ), dtype = np.complex128 ) 
    #filter and fill the mesh
    A = spectr #i am using the same letters of the formula
    r = freq * T 
    
    for i in range(0,N):
        r_zero =  interest_freq[i] * T 
        for j in range(0,M):
            rdot = range_rdot[j]
            m =  np.int( 2 * rdot)#np.abs( np.rint( ( 2 * rdot ) ) )#maximum_m_idx
            if( np.abs( rdot ) > 0.01 ):

                factor = 1 / scipy.sqrt( 2 * rdot )

                kmin = center - int( (m / 2) )
                kmax = center + int( (m / 2) )

                element = 0
                for k in range(kmin, kmax +1): #( kmin -1, kmax  )
                    q_k =  r_zero - np.int( r[k] ) 

                    exponential = scipy.exp(1j * np.pi * q_k**2 /  rdot  )

                    Yk = scipy.sqrt( 2 * rdot ) * q_k
                    Zk = scipy.sqrt( 2 * rdot ) * ( q_k + rdot )
                    [SZ, CZ] = special.fresnel(Zk)
                    [SY, CY] = special.fresnel(Yk)
                    weight = SZ - SY - 1j * (CY - CZ)

                    element = element + A[k] * weight * exponential * factor 

                mesh[i,j] =  element
        print (  str(round( 100 * i / N, 2 ) ) + "%")
    rdot_center = np.argmin(np.abs(range_rdot))    
    np.transpose(mesh)[rdot_center] = interest_spectr

           
    if vi == True :
        accel_utils.plot(interest_freq, range_rdot, np.abs( mesh ) ) 
    
    
    maximum = np.unravel_index(np.argmax(abs(mesh), axis=None), mesh.shape)
    # print( "max index on the matrix: " + str( maximum ))
    
    return interest_freq[ maximum[0] ] , range_rdot[maximum[1]] / T**2 #+ int( delta_rdot**-1 ) 




#this function get the power, find an interesting point to start to seach, cut the FFT and
#call call the accelsearch method to find f and fdot
def find_frequency_and_fdot( times, signal, vi = False, interbin = False, delta_r = 1, delta_rdot = 1 ):

    N, dt, T, spectr, freq = accel_utils.get_info( times, signal)
    uncert_f = accel_utils.get_uncert_f( spectr, freq ) #and uncertain deltaR ?
    # print(uncert_f)
    if interbin:
        freq, spectr = accel_utils.my_interbin_fft(freq, spectr)
        N = 2 * N
        #T = 2 * T
    # if delta_r < 1:
    #     T, N, freq, spectr = accel_utils.fft_interpolation(spectr, freq, T, dt, delta_r, uncert_f )  
    int_spectr, int_freq  = accel_utils.interest_search( spectr, freq, uncert_f )
    #print(int_spectr.size, int_freq.size)
    #print(int_freq[int(len(int_spectr) / 2)])
     
    return accelsearch( int_spectr, int_freq, T, vi , delta_rdot)

def plot_f_fdot( times, signal, interbin = False, delta_r = 1, delta_rdot = 1 ):
    flag = interbin
    r = delta_r
    rdot = delta_rdot
    return find_frequency_and_fdot(times, signal, vi = True, interbin = flag, delta_r = r, delta_rdot = rdot )
    
    
    
    ###############################################################################################################################
    
    import numpy as np
import scipy
from scipy import special
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from mpl_toolkits.mplot3d import axes3d
from scipy import signal
from scipy.fftpack import fft, fftshift
from scipy.fftpack import fft, fftshift
import sys
import accel_utils


def nearest_value_idx( array, value ):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


def nearest_freq( power , freq, pow_value ):
    idx = nearest_value_idx( power, pow_value )
    return freq[idx]


def get_power( spectr ):
    return np.abs( spectr )**2


def get_norm_power( spectr ): 
    Nph = np.abs( spectr[0] )
    return 2 * get_power ( spectr ) / Nph


def get_dt( times ):
    return times[1] - times[0]

# it is focused on the postive frequencys
def get_uncert_f( spectr, frequency ):
    N = len(frequency)
    power = get_norm_power( spectr )
    _, bins = np.histogram( power , bins = 50 )
    return nearest_freq( power[0:int(N / 2)], frequency[0:int(N / 2)], np.median( bins ) )


#it is usefull to study only a little portion of the f-fdot plane
def interest_search(  spectr, freq, uncert_f ):
    # i should add a delta_R code to choose the best min freq and max freq
    N = len( freq )
    delta =  5 * np.int( N / ( 10 ** (np.log10( N ) -2 ) ) + 1) 
    center = nearest_value_idx( freq, uncert_f ) #int(N / 2) #
    lside = center - delta
    rside = center + delta + 1
    return spectr[ lside : rside ], freq[ lside : rside ] 

def get_info( times, signal ):
    N = len( times )
    dt = get_dt( times )
    T = N * dt
    spectr = np.fft.fft( signal )[0 : int( N / 2) ] 
    freq = np.fft.fftfreq( N , dt)[0 : int( N / 2) ] 
    return N, dt, T, spectr, freq



def plot( freq, range_rdot, mesh):
    # X = freq
    # Y = range_rdot
    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    # X, Y = np.meshgrid( X, Y )

    # # Plot a basic wireframe.
    # ax.plot_wireframe(Y, X,  np.transpose( mesh ), rstride=10, cstride=10)
    # plt.show()

    #2d projection
    plt.imshow( np.transpose( abs(mesh) ), aspect = 'equal'  )
    plt.show()



def interbin_fft(freq, fft):
    """
    Examples
    --------
    >>> freq = [-1, -0.5, 0, 0.5, 1]
    >>> fft = np.array([1, 0, 1, 0, 1])
    >>> f, F = interbin_fft(freq, fft)
    >>> np.allclose(f, [-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1])
    True
    >>> pi_4 = np.pi / 4
    >>> np.allclose(F, [1, -pi_4, 0, pi_4, 1, -pi_4, 0, pi_4, 1])
    True
    """
    freq = np.array(freq)
    fft = np.array(fft)
    order = np.argsort(freq)

    freq = freq[order]
    fft = fft[order]

    new_freqs = np.linspace(freq[0], freq[-1], 2 * len(freq) - 1)
    new_fft = np.zeros_like(new_freqs, dtype=type(fft[0]))
    new_fft[::2] = fft
    new_fft[1::2] = (fft[1:] - fft[:-1]) * np.pi / 4 
    return new_freqs, new_fft


