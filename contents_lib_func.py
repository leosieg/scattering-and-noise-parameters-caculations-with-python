lib_conN.py
    #
    def Nnpc_xy -> port connection x-port -> y-port, signal and noise
    #
    def Nnpc_x -> port connection at m-port, signal and noise
    #
lib_conS.py
    #
    def Snpc_xy -> port connection x-port -> y-port, signal only
    #
    def Snpc_x -> port connection at m-port, signal only
    #
lib_line.py
    #    
    def line_s -> mismatched homogeneous series line
    #
    def line_p -> parallel stub line
    #
    def gamma_z -> reflection coefficient for z
    #
    def gamma_y -> reflection coefficient for y
    #
    def y_parallel -> parallel admittance
    #
    def z_series -> series impedance
    #
    def jomegalc -> factor for normalizing inductance or capacitance
    #
    def p3sp -> three port splitter, common earth
    #
    def wilkinson_s -> 50 Ohm Wilkinson splitter 
    #
lib_Ncircuit
    #
    def lna_circuit_N -> circuit structure low noise amplifier 
    #                    for signal and noise analysis
    #
    def amp_circuit_N -> circuit structure balanced amplifier
    #                    for signal and noise analysis
    #
    def amp_circuit_N02 -> circuit structure double balanced amplifier
    #                      for signal and noise analysis
    #
    def cascode_circuit_N -> circuit structure cascode circuit
    #                        for signal and noise analysis 
    #
lib_noise
    #
    def ncir -> noise circles plot function
    #
    def ncirp -> noise circles and paraboloid plot function
    #    
lib_Scircuit
#
    def lna_circuit_S -> circuit structure low noise amplifier 
    #                    for signal analysis only
    #
    def amp_circuit_S -> circuit structure balanced amplifier
    #                    for signal analysis only
    #
    def amp_circuit_S02 -> circuit structure double balanced amplifier
    #                      for signal analysis only
    #
    def cascode_circuit_N -> circuit structure cascode circuit
    #                        for signal analysis only 
    #
    def filter_circuit_S -> circuit structure filter circuit
    #                       for signal analysis only 
    #
lib_smith_stab
    #
    def smith_z -> Smith chart impedance(z) version
    #
    def smith_y -> Smith chart admittance(y) version
    #
    def ssks -> stability circle source plane
    #
    def sskl -> stability circle load plane 
    #
    def smith3_z -> Smith chart impedance(z) version for 3d projection
    #
    def smith3_y -> Smith chart admittance(y) version for 3d projection   
    #
lib_spline
    #
    def s2spl_S -> spline approximation for data in tables *.s2p
    #              signal only
    #
    def s2spl_N -> spline approximation for data in tables *.s2p
    #              signal and noise
    #
lib_t2tot3
    #
    def t2tot3ce_N -> source free two port -> three port, common earth
    #                 signal and noise
    #
    def t2tot3ce_S -> source free two port -> three port, common earth
    #                 signal only
    #   
    def t2tot3ef_N -> source free two port -> three port, earth free
    #                 signal and noise
    #
    def t2tot3ef_S -> source free two port -> three port, earth free
    #                 signal only
    #  