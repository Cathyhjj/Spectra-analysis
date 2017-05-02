
# coding: utf-8
#RIXS_Data_Analysis

"""
Exprimental RIXS Data Anaylysis, ESRF ID26 
--------- version(1.0) --------- 
HAVE FUN WITH YOUR DATA ANALYSIS!
02-05-2017 -- Juanjuan Huang(Cathy). 

This includes: 
(I)  XANES data processing
    1.1 Averaged/summed XANES plotting with interpolation for incident energy
    1.2 XANES area normalization (normalized to whole area or specified tail region)
(II) RIXS data processing
    2.1 2D/3D RIXS plane plotting with interpolation for both incident energy and emission energy
        2.1.1 concentration correction
        2.1.2 IE versus EE plotting 
        2.1.3 IE versus ET plotting 
        2.1.4 RIXS plane CIE CEE CET cuts and integration plotting
    2.2 Averaging and merging for RIXS planes
(III)Save the data into RIXS txt files so that can be further used by other software, e.g, Matlab

-----
To do
[] Center_of_Mass of the peaks
[] Peaks fitting
[] Save into a txt file
"""

import numpy as np
import matplotlib.pyplot as plt
from silx.io.specfile import SpecFile
import scipy.interpolate as interp
import scipy.ndimage as nd
from scipy import signal
import os

class DataAnalysis(object):
    '''
 |   class DataAnalysis(object)
 |
 |   This includes: 
 |   1. Averaged/summed XANES plotting with interpolation for incident energy
 |   2. 2D/3D RIXS plane plotting with interpolation for both incident energy and emission energy
 |      2.1 concentration correction
 |      2.2 IE versus EE plotting 
 |      2.3 IE versus ET plotting 
 |      2.4 RIXS plane CIE CEE CET cuts and integration plotting
 |   3. Averaging and merging for RIXS planes
 |   4. Save the data into RIXS txt files so that can be further used by other software, e.g, Matlab
 |
 |
 |  Methods
 |  ----------
 |  (for the __new__ method; see Notes below)
 |
 |  XANES_data(): get XANES merged data ndarray from SPEC file
 |      return 1d ndarray [incident energy, intensity]
 |
 |  RIXS_data() : To get RIXS data ndarray from SPEC file
 |      return data ndarray [incident energy, emission energy, intensity]
 |
 |  RIXS_merge: To merge (sum up/average different RIXS data ndarray)
 |      return data ndarray [incident energy, emission energy, intensity]
 |
 |  RIXS_display : To plot RIXS planes
 |      return None
 |
 |  RIXS_cut : To do CIE, CET, CEE cuts
 |      return CIE, CET, CEE data ndarray [incident energy/energy transfer, intensity]
 |
 |  RIXS_integration : Integration along incident energy and energy transfer
 |      return integrated data ndarray [incident energy/energy transfer, intensity]
 |
 |  Parameters
 |  ----------
 |  path : the filepath of Specfile
    '''
    
    def __init__(self, path):
        self.path = path
        self.sf = SpecFile(path)
    
    def XANES_data(self, firstScan, lastScan, interp_npt_1eV = 20, method = 'average', savetxt = False):
        """
        To get XANES merged data ndarray from SPEC file
        The incident energy for scans can be different

        Parameters
        ----------
        firstScan : the index of first scan, e.g, 71 corresponding to fscan '72.1'
        lastScan : the index of first scan
        interp_npt_1eV : the number of interpolation points for 1 eV, default: 20 points for 1 eV
                         e.g, Incident Energy: 6535 eV - 6545 eV, 11 eV, 115 points, -----> 220 points
        method : 'average' or 'sum' for intensity
        savetxt: default True, save the ET, EE data as folders 
        Returns
        -------
        out : A 1d data ndarray [incident_Energy_interp, XANES_merge_inten]
              incident_Energy_interp -----> interpolated incident energy
              XANES_merge_inten -----> interpolated intensity
    """
        # Each scan has different incident energy points
        # this step finds the highest incident energy corresponding scan
        #             and the lowest incident energy corresponding scan
        energy_checkmin_list = []
        energy_checkmax_list = []
        for n in range(firstScan, lastScan + 1): 
            IE_min_check = self.sf[n].data_column_by_name('arr_hdh_ene')[0]        
            IE_max_check = self.sf[n].data_column_by_name('arr_hdh_ene')[-1]
            energy_checkmin_list.append(IE_min_check)
            energy_checkmax_list.append(IE_max_check)
        energy_min_index = energy_checkmin_list.index(min(energy_checkmin_list)) + firstScan
        energy_max_index = energy_checkmax_list.index(max(energy_checkmax_list)) + firstScan

        # Find the energy span of incident energy (For later interpolation)
        # e.g, incident energy range : 4.987654 KeV - 4.987987 KeV
        # will be approximated as 4.9878 KeV for minimum incident Energy 
        #                         4.9879 KeV for maxmum incident Energy
        # I do in such a way to ensure the interpolation points lie always inside the experimental incident energy range
        # 4.9878 > 4.987654 while 4.9879 < 4.987987
        incident_Energy_min = round(self.sf[energy_min_index].data_column_by_name('arr_hdh_ene')[0]*10000+1)/10000
        incident_Energy_max = round(self.sf[energy_max_index].data_column_by_name('arr_hdh_ene')[-1]*10000-1)/10000

        # Find the energy span of incident energy
        incident_Energy_Span = incident_Energy_max - incident_Energy_min

        # Define the total points of interpolation for incident energy
        # default: 20 points for 1 eV
        incident_Energy_interp_npt = int(round(incident_Energy_Span*1000) * interp_npt_1eV)

        # Fisrt do the incident energy 1d interpolation
        # Creat empty arrays filled with NaN for different XANES scans 
        # which has the shape (scan total numbers, incident_Energy_interp_npt)
        XANES_inten_array = np.zeros(((lastScan - firstScan + 1,incident_Energy_interp_npt))) 
        XANES_inten_array[:] = np.nan

        # We then fill the empty array with interpolated concentration corrected intensity
        # Use scipy.interpolate.interp1d to interpolate the points
        for n in range(firstScan, lastScan + 1):
            incident_Energy = self.sf[n].data_column_by_name('arr_hdh_ene')
            # corresponding intensity
            inten = self.sf[n].data_column_by_name('det_dtc')/self.sf[n].data_column_by_name('I02') # Normalized to I02
            # Define our interp1d function for incident energy, fill the outbound value with 0
            f_interp = interp.interp1d(incident_Energy, inten, bounds_error = False)
            # Find our interpolated incident energy and corresponding intensity
            incident_Energy_interp = np.linspace(incident_Energy_min, incident_Energy_max, incident_Energy_interp_npt)
            inten_interp = f_interp (incident_Energy_interp)
            # Fill the empty array with the interpolated data
            # Now we only interpolate in incident axis, not in emission axis
            XANES_inten_array[n-firstScan,:] = inten_interp

        if method == 'average':
            # To average all the intensities for different scans, we ignore the nan data
            XANES_merge_inten = np.nanmean(XANES_inten_array, axis = 0)
        elif method == 'sum':
            # To average all the intensities for different scans, we ignore the nan data
            XANES_merge_inten = np.nansum(XANES_inten_array, axis = 0)

        # Put incident energy and merged intensity into XANES data array
        dataArray_XANES = np.array([incident_Energy_interp, XANES_merge_inten]) 

        if savetxt == True:
            pass
        
        return dataArray_XANES
    

    def XANES_normalize(self, XANES_data, normalized_starting_energy = None):
        """
        To do normalization for XANES
        Parameters
        ----------
        XANES_data : the XANES_data output ndarray
        normalized_starting_energy: An energy number, e.g, Incident Energy: 6600 eV, 
                                              -----> will do normalization from 6600 eV to the end of tail feature
                                    Not difine-----> Normalization to the whole area
        Returns
        -------
        out : A data ndarray [energy, normlized_intensity]

        """
        if normalized_starting_energy == None:
            normalized_starting_energy = XANES_data[0][0] * 1000
        postedge_index = np.where(XANES_data[0] >= normalized_starting_energy/1000)
        postedge_intensity = XANES_data[1][postedge_index[0]]
        # Calculate the tail edge area
        tail_edge_area = np.trapz(postedge_intensity, dx=1)
        # Normalization to the whole area
        norm_intensity = XANES_data[1]/tail_edge_area
        norm_dataArray = np.array([XANES_data[0],norm_intensity])
        return norm_dataArray
    

    def RIXS_data(self,firstScan, lastScan, concCorrecScan, interp_npt_1eV = 20, choice = 'EE', savetxt = False):
        """
        To get RIXS data ndarray from SPEC file

        Parameters
        ----------
        firstScan : the index of first scan, e.g, 71 corresponding to fscan '72.1'
        lastScan : the index of first scan
        concCorrecScan : the index of concentration correction scan, normally the one after last RIXS scan
        interp_npt_1eV : the number of interpolation points for 1 eV, default: 20 points for 1 eV
                         e.g, Incident Energy: 6535 eV - 6545 eV, 11 eV, 115 points, -----> 220 points
                              Emitted  Energy: 5890 eV - 5905 eV, 15 eV, 76  points, -----> 300 points
        choice : 'EE': get -----> incident energy & emission energy plotting
                 'ET': get -----> energy transfer & emission energy plotting
        savetxt: default True, save the ET, EE data as folders 
        Returns
        -------
        if choice = 'EE'
        out : ndarray, A data list [EE_XX, EE_YY, EE_MDfci_correc_inten_2dinterp]
              EE_XX -----> interpolated incident energy ndarray
              EE_YY -----> interpolated emission energy ndarray
              EE_MDfci_correc_inten_2dinterp -----> interpolated intensity ndarray

        if choice = 'ET'
        out : ndarray, A data list [ET_XX, ET_YY, ET_MDfci_correc_inten_2dinterp]
              ET_XX -----> interpolated incident energy ndarray
              ET_YY -----> interpolated energy transfer ndarray
              ET_MDfci_correc_inten_2dinterp -----> interpolated intensity ndarray
    """

        # Extract emission energy from SPEC file
        emission_Energy = np.array([self.sf[i].data_column_by_name('xes_en')[1] for i in range(firstScan,(lastScan+1))])    

        # Find the energy span of incident energy
        incident_Energy_min = round(self.sf[firstScan].data_column_by_name('arr_hdh_ene')[0]*10000+1)/10000
        incident_Energy_max = round(self.sf[firstScan].data_column_by_name('arr_hdh_ene')[-1]*10000-1)/10000

        # and emitted energy
        emission_Energy_min = round(self.sf[firstScan].data_column_by_name('xes_en')[0]*10000+1)/10000
        emission_Energy_max = round(self.sf[lastScan].data_column_by_name('xes_en')[0]*10000-1)/10000

        # Find the energy span of incident energy
        incident_Energy_Span = incident_Energy_max - incident_Energy_min
        # And emitted energy
        emission_Energy_Span = emission_Energy_max - emission_Energy_min

        # Define the total points of interpolation for incident energy
        incident_Energy_interp_npt = int(round(incident_Energy_Span*1000)*interp_npt_1eV)
        # And emission energy
        emission_Energy_interp_npt = int(round(emission_Energy_Span*1000)*interp_npt_1eV)

        # Collect concentration correction intensity into an array
        concCorrec_inten = self.sf[concCorrecScan].data_column_by_name('det_dtc')/self.sf[concCorrecScan].data_column_by_name('I02') # Normalized to I02

        # Fisrt do the incident energy 1d interpolation
        # Creat empty arrays filled with 0 for RIXS intensity 
        # which has the shape (emission Energy (scan total numbers), incident_Energy_interp_npt)
        MDfci_correc_inten = np.zeros(((len(emission_Energy),incident_Energy_interp_npt)))        
        # We then fill the empty array with interpolated concentration corrected intensity
        # Use scipy.interpolate.interp1d to interpolate the points
        for n in range(firstScan, lastScan + 1):
            incident_Energy = self.sf[n].data_column_by_name('arr_hdh_ene')
            # To do concentration correction for intensity
            correc_inten = self.sf[n].data_column_by_name('det_dtc')/(self.sf[n].data_column_by_name('I02')*concCorrec_inten[n-firstScan]) #Normalized to I02
            # Define our interp1d function for incident energy, fill the outbound value with 0
            f_interp = interp.interp1d(incident_Energy, correc_inten, bounds_error = False, fill_value = 0)
            # Find our interpolated incident energy and corresponding intensity
            incident_Energy_interp = np.linspace(incident_Energy_min, incident_Energy_max, incident_Energy_interp_npt)
            correc_inten_interp = f_interp (incident_Energy_interp)
            # Fill the empty array with the interpolated data
            # Now we only interpolate in incident axis, not in emission axis
            MDfci_correc_inten[n-firstScan,:] = correc_inten_interp

        # After doing 1D interpolation for incident energy
        # Now we are going to do 2D interpolation for emission energy, using scipy.interpolate.interp2d
        # Define our interp2d function
        f_interp2d = interp.interp2d(incident_Energy_interp, emission_Energy, MDfci_correc_inten)
        # Define our new emission energy
        emission_Energy_interp = np.linspace(emission_Energy_min, emission_Energy_max, emission_Energy_interp_npt)
        # Get interpolated new intensity array (emission energy interpolated)
        EE_MDfci_correc_inten_2dinterp = f_interp2d(incident_Energy_interp, emission_Energy_interp)

        # Define Grids, EE_XX: incident energy array, EE_YY: emission energy array
        EE_XX, EE_YY = np.meshgrid(incident_Energy_interp, emission_Energy_interp)

        # Put all the data into a list
        dataArray_EE = np.array([EE_XX, EE_YY, EE_MDfci_correc_inten_2dinterp])

        # -------------- RIXS Energy Transfer - Incident Energy plotting 
        # When it comes to ET, the length of new y axis(energy transfer) change
        energy_transfer_min = min(EE_XX[0,:])-max(EE_YY[:,0])
        energy_transfer_max = max(EE_XX[0,:])-min(EE_YY[:,0])  
        energy_transfer_length = sum(EE_XX.shape)-1
        # Define our energy transfer axis
        energy_transfer = np.linspace(energy_transfer_min,energy_transfer_max,energy_transfer_length)

        # Define our new intensity empty array: ET_MDfci_correc_inten_2dinterp
        # And the new empty array has a shape of (sum(EE_XX.shape)-1, EE_XX.shape[1])
        ET_MDfci_correc_inten_2dinterp = np.zeros((energy_transfer_length, EE_XX.shape[1]))
        # Make all the numbers into NaN
        ET_MDfci_correc_inten_2dinterp[:] = np.nan
        offset = (EE_YY.shape[0]-1,0)

    #     #---------------- AFFINE TRANSFORM METHOD ----------------
    #     # We want new_y = x-y, new_x = x
    #     # Thus we need the transform matrix : [-1,1]
    #     #                                     [0, 1]
    #     transform_matrix = np.array([[-1,1],[0,1]])
    #     ET_MDfci_correc_inten_2dinterp = nd.interpolation.affine_transform(EE_MDfci_correc_inten_2dinterp, transform_matrix, 
    #                                                               output_shape = ET_MDfci_correc_inten_2dinterp.shape, offset = offset)

        # Define ET Grids, ET_XX: incident energy array, ET_YY: energy transfer array
        ET_XX, ET_YY = np.meshgrid(incident_Energy_interp, energy_transfer)

        #---------------- LOOP METHOD ----------------
        for j in range(EE_XX.shape[0]):
            for i in range(EE_XX.shape[1]):
                # Fill the new empty array with the EE intensity
                ET_MDfci_correc_inten_2dinterp[i-j+offset[0],i] = EE_MDfci_correc_inten_2dinterp[j,i]


        # Put all the data into a list
        dataArray_ET = np.array([ET_XX, ET_YY, ET_MDfci_correc_inten_2dinterp])

        if savetxt == True:
            # Save file: Creat EE and ET folders in the compound file folder
            EE_path = fileFolder + 'Emission_Energy'
            ET_path = fileFolder + 'Energy_Transfer'
            if not os.path.exists(EE_path):
                os.makedirs(EE_path)
                np.savetxt(EE_path + '_EE.txt', dataArray_EE, fmt = '%.10f')
            if not os.path.exists(ET_path):
                os.makedirs(ET_path)
                np.savetxt(EE_path + '_ET.txt', dataArray_ET, fmt = '%.10f')

        if choice == 'EE':
            return dataArray_EE
        if choice == 'ET':      
            return dataArray_ET

    def RIXS_merge(self, scansets, choice = 'sum'):
        """
        To merge (sum up/average different RIXS data ndarray)

        Parameters
        ----------
        scansets : put all the RIXS_data output file that want to mergy into a list
                   e.g. [dataArray1, dataArray2, dataArray3]
        choice : 'sum':     get -----> summed intensity
                 'average': get -----> averaged intensity
        Returns
        -------
        if choice = 'sum':
        out : A new ndarray dataArray with summed intensity 
              -----> [averaged_XX,averaged_YY,summed_intensity]

        if choice = 'average':
        out : A new ndarray dataArray with summed intensity
              -----> [averaged_XX,averaged_YY,averaged_intensity]
        """
        # scansets = np.array(scansets)
        summed_XX = 0
        summed_YY = 0
        summed_intensity = 0
        # Sum up incident energy array, emission energy array and intensity array separately
        for nscan in scansets:
            summed_XX += nscan[0]
            summed_YY += nscan[1]
            summed_intensity += nscan[2]

        # To do the average = sum/scansets
        averaged_XX = summed_XX/len(scansets)
        averaged_YY = summed_YY/len(scansets)
        averaged_intensity =summed_intensity/len(scansets)

        # Reput the data into a dataArray
        summed_dataArray = np.array([averaged_XX,averaged_YY,summed_intensity])
        averaged_dataArray = np.array([averaged_XX,averaged_YY,averaged_intensity])

        if choice == 'sum':
            return summed_dataArray
        if choice == 'average':
            return averaged_dataArray
    
    def RIXS_display(self, dataArray, title = 'RIXS', choice = 'EE', mode = '2d'):
        """
        To plot RIXS planes

        Parameters
        ----------
        dataArray: the return ndarray [XX, YY, intensity] from RIXS_data
        title: plotting title, str
        x_lim: set limit for x axis, e.g, (6537,6544) stands for energy range from 6537-6544 eV
        x_lim: set limit for y axis, e.g, (5892,5902)
        choice: 'EE' emission energy(default) or 'ET' energy transfer
        mode: '2d' 2D plotting(default), '3d' 3D plotting
        Returns
        -------
        out : ndarray
            Array of zeros with the given shape, dtype, and order.
    """
        if mode == '2d':
            # ------------ CONTOURF METHOD ------------
            # Transfer the energy from KeV scale into eV scale
            XX = dataArray[0] * 1000
            YY = dataArray[1] * 1000
            intensity = dataArray[2]
            plt.figure()
            plt.contourf(XX, YY, intensity, 20)
            plt.colorbar()
            plt.contour(XX, YY, intensity, 20, linewidths = 0.5, colors='black')
            plt.title(title)
            #    plt.axis('equal')

            # ------------ IMSHOW METHOD ------------       
            #         plt.imshow(scan1_ET[2],origin='lower', extent=[scan1_ET[0][0,0]*1000, scan1_ET[0][0,-1]*1000, scan1_ET[1][0,0]*1000, scan1_ET[1][-1,0]*1000], interpolation="nearest")
            #         plt.axis('equal')
            #         plt.figsize = (scan1_ET[2].shape)
            #         plt.show()

            plt.xlabel('Incident Energy [eV]')
            if choice == 'EE':
                plt.ylabel('Emitted Energy [eV]')
                # plt.ylim(y_lim)
                # plt.xlim(x_lim)
                # x_lim = (6537,6544),y_lim=(5892,5902), 
            elif choice == 'ET':
                plt.axis('equal')
                plt.ylabel('Energy Transfer [eV]')
        elif mode == '3d':
            from mpl_toolkits.mplot3d import Axes3D
            from matplotlib import cm
            from matplotlib.ticker import LinearLocator, FormatStrFormatter
            XX = dataArray[0] * 1000
            YY = dataArray[1] * 1000
            intensity = dataArray[2]
            fig = plt.figure()
            ax = fig.gca(projection='3d')
            plotting_3d = ax.plot_surface(XX, YY, intensity, cmap=cm.coolwarm, linewidth=0, antialiased=False)
            # Customize the z axis.
            ax.zaxis.set_major_locator(LinearLocator(10))
            ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
            # Add colorbar
            fig.colorbar(plotting_3d, shrink=0.5, aspect=5)
        return plt.show()
    
    def RIXS_cut(self, dataArray, choice, cut):
        """
        To do CIE, CET, CEE cuts(choice, cut)
        NOTICE: Choose ET dataArray for CIE & CET
                Choose EE dataArray for CEE

        Parameters
        ----------
        dataArray: the RIXS_data output file
        choice: 'CIE'-- Constant incident energy cut
                'CET'-- Constant energy transfer cut
                'CEE'-- Constant emission energy cut
        cut: the energy (eV) you want to cut, e.g., 6530 eV

        Returns
        -------
        out : 
        ndarray
 |      CIE, CET, CEE data ndarray [incident energy/energy transfer, intensity]
    """

        # Convert all the NaN to numbers (the 2d interpolation will raise errors when existing NaN)
        new_inten = np.nan_to_num(dataArray[2])
        # 2d interpolation
        cut_interp2d = interp.interp2d(dataArray[0][0,:], dataArray[1][:,0], new_inten)
        # Convert KeV into eV
        cut = cut/1000
        if choice == 'CIE':
            # Find the interpolated intensity
            cut_intensity = cut_interp2d(cut, dataArray[1][:,0])
            # Plotting
            plt.plot(dataArray[1][:,0]*1000, cut_intensity)
            plt.title('CIE')
            plt.xlabel('Energy transfer')
            plt.ylabel('Arbitrary Intensity')
            plt.show()
            CIE_dataArray = np.array([dataArray[1][:,0]*1000, cut_intensity[0]])
            return CIE_dataArray
        elif choice == 'CET':
            # Find the interpolated intensity
            cut_intensity = cut_interp2d(dataArray[0][0,:],cut)
            # Plotting
            plt.plot(dataArray[0][0,:]*1000,cut_intensity)
            plt.title('CET')
            plt.xlabel('Incident Energy')
            plt.ylabel('Arbitrary Intensity')
            plt.show()
            CET_dataArray = np.array([dataArray[0][0,:]*1000,cut_intensity])
            return CET_dataArray
        elif choice == 'CEE':
            # Find the interpolated intensity
            cut_intensity = cut_interp2d(dataArray[0][0,:],cut)
            # Plotting
            plt.plot(dataArray[0][0,:]*1000,cut_intensity)
            plt.title('CEE')
            plt.xlabel('Incident Energy')
            plt.ylabel('Arbitrary Intensity')
            plt.show()
            CEE_dataArray = np.array([dataArray[0][0,:]*1000,cut_intensity])
            return CEE_dataArray
    
    def RIXS_integration(self, dataArray):
        """
        Integration along incident energy and energy transfer

        Parameters
        ----------
        dataArray: the RIXS_data return data ndarray

        Returns
        -------
        out :     
        ndarray
        integrated data ndarray [incident energy/energy transfer, intensity]


    """
        # integration for incident energy ---> Conventional XANES
        sumIntensity_IE = np.nansum(dataArray[2],axis = 0)
        # integration for incident energy ---> Conventional XANES
        sumIntensity_ET = np.nansum(dataArray[2],axis = 1)
        # plot both figures with same intensity scale
        fig, ax = plt.subplots(nrows=1, ncols=2,figsize=(12,4),sharey = 'all')
        ax[0].plot(dataArray[0][0,:]*1000,sumIntensity_IE)
        ax[0].set_xlabel('Incident Energy [eV]')
        ax[0].set_ylabel('Integrated intensity')
        ax[1].plot(dataArray[1][:,0]*1000,sumIntensity_ET)
        ax[1].set_xlabel('Energy Transfer [eV]')
        ax[1].set_ylabel('Integrated intensity')
        plt.setp(ax[1].get_yticklabels(), visible = True)
        plt.show()
        integration_dataArray = np.array([[dataArray[0][0,:]*1000,sumIntensity_IE],[dataArray[1][:,0]*1000,sumIntensity_ET]])
        return integration_dataArray
