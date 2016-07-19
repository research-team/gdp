import nest
import numpy as np
import pylab as pl
import nest.raster_plot

nest.ResetKernel()

# neuron_model = "iaf_psc_alpha"
# neuron_model = "iaf_psc_delta"
# neuron_model = "iaf_psc_exp"
# neuron_model = "iaf_cond_exp"
# neuron_model = "iaf_cond_alpha"
# neuron_model = "aeif_cond_alpha"
# neuron_model = "mat2_psc_exp"
# neuron_model = "hh_psc_alpha"
neuron_model = "hh_cond_exp_traub"

# spike detector
sd = nest.Create("spike_detector")
nest.SetStatus(sd, {"label": "spikes", "withtime": True, "withgid": True, "to_file": False})

# display recordables for illustration
print neuron_model, ' recordables: ', nest.GetDefaults(neuron_model)['recordables']

# create neuron and multimeter
'''
V_m double - Membrane potential in mV
V_T double - Voltage offset that controls dynamics. For default parameters, V_T = -63mV results in a threshold around -50mV.
E_L double - Resting membrane potential in mV.
C_m double - Capacity of the membrane in pF
g_L double - Leak conductance in nS.
tau_m double - Membrane time constant in ms.
t_ref double - Duration of refractory period in ms.
V_th double - Spike threshold in mV.
V_reset double - Reset potential of the membrane in mV.
tau_syn_ex double - Rise time of the excitatory synaptic alpha function in ms.
tau_syn_in double - Rise time of the inhibitory synaptic alpha function in ms.
I_e double - Constant external input current in pA.
V_min double - Absolute lower value for the membrane potential.

E_ex double - Excitatory synaptic reversal potential in mV.
E_in double - Inhibitory synaptic reversal potential in mV.
E_Na double - Sodium reversal potential in mV.
g_Na double - Sodium peak conductance in nS.
E_K double - Potassium reversal potential in mV.
g_K double - Potassium peak conductance in nS.
'''

if neuron_model == "hh_cond_exp_traub":
    n = nest.Create(neuron_model, params={'V_m': -70.0, 'E_L': -70.0, 'V_T': -63.0, 'g_L': 100.0, 'C_m': 30.0,
                                          'tau_syn_ex': 5.0, 'tau_syn_in': 15.0, 'I_e': 0.0,
                                          'E_ex': 0.0,
                                          'E_in': -80.0,
                                          'E_Na': 50.0,
                                          'g_Na': 20000.0,
                                          'E_K': -90.0,
                                          'g_K': 6000.0
                                          })
else:
    n = nest.Create(neuron_model, params={'tau_syn_ex': 1.0})

m = nest.Create('multimeter', params={'withtime': True, 'interval': 0.1, 'record_from': ['V_m', 'I_Na', 'I_K', 'I_Cl']})

# Create spike generators and connect
gex = nest.Create('spike_generator', params={'spike_times': np.array([1.0])})
gin = nest.Create('spike_generator', params={'spike_times': np.array([6.0])})
# gex = nest.Create('spike_generator', params={'spike_times': np.array([5.0, 10, 15])})
# gin = nest.Create('spike_generator', params={'spike_times': np.array([15.0, 25.0, 55.0])})

nest.Connect(gex, n, syn_spec={'weight': 40.0})  # excitatory
nest.Connect(gin, n, syn_spec={'weight': -40.0})  # inhibitory
nest.Connect(m, n)
nest.Connect(n, sd)

# simulate
nest.Simulate(10)

# obtain and display data
events = nest.GetStatus(m)[0]['events']
t = events['times']

pl.subplot(211)
pl.plot(t, events['V_m'])
pl.ylabel('Membrane potential [mV]')

pl.subplot(212)
pl.plot(t, events['I_Na'], t, events['I_K'], t, events['I_Cl'])
pl.ylabel('Currents [nA]')
pl.legend(('I_Na', 'I_K', 'I_Cl'))

nest.raster_plot.from_device(sd)
nest.raster_plot.show()
