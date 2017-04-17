import DataAnalysis
import matplotlib.pyplot as plt


# XANES test
compound7 = DataAnalysis.DataAnalysis('/Users/brookcathy/Desktop/juhuang/compound7/Compound_7')
energy1, intensity1 = compound7.XANES_data(7,10, method = 'sum')
energy2, intensity2 = compound7.XANES_data(3,5, method = 'sum')
energy3, intensity3 = compound7.XANES_data(3,62)
fig = plt.figure(1)
plot1 = plt.plot(energy1,intensity1,energy2,intensity2,energy3,intensity3)
plt.show()

# RIXS test
scan1 = compound7.RIXS_data(71, 146, 147)
compound7.RIXS_display(scan1)
scan1_ET = compound7.RIXS_data(71, 146, 147, choice = 'ET')
compound7.RIXS_display(scan1_ET,choice='ET')

# RIXS Cut test
compound7.RIXS_cut(scan1_ET,'CIE',6540)
compound7.RIXS_cut(scan1_ET,'CET',645)
compound7.RIXS_cut(scan1_ET,'CET',650)
compound7.RIXS_cut(scan1,'CEE',5898)
# Integration Test
compound7.RIXS_integration(scan1_ET)
