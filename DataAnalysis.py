# coding: utf-8
#RIXS_Data_Analysis

"""
Exprimental RIXS Data Anaylysis, ESRF ID26 
--------- version(1.1) ------------ 
HAVE FUN WITH YOUR DATA ANALYSIS!
      Juanjuan Huang(Cathy) 
-----------------------------------
"""

# LOGBOOK
# 20170622 -- update : Choice of concen.correc in RIXS planes, RIXS_data_constant_ET() method, Normalization plotting 
# 20170615 -- update : RIXS_normalization() method
# 20170613 -- update : skip problematic scans
# 20170604 -- update : RIXS_imshow() additional method
# 20170528 -- update : Radiation damage special average, Radiation_damage() method
# 20170506 -- update : Add general functions: saveFile(), normalize_toArea()
# 20170503 -- update : Adding XANES_find_peaks() method
# 20170502 -- update : Adding Normalization of XANES data
# 20170419 -- First version

# To do
# 1. Checker
#    2.1 Check whether intensity is wrong
#    2.2 Check concentration correction

# 2. XANES plotter, XANES KeV, eV choice and change other XANES method units to eV

import numpy as np
import matplotlib.pyplot as plt
from silx.io.specfile import SpecFile
import scipy.interpolate as interp
import scipy.ndimage as nd
import os
from scipy import signal
from matplotlib import cm

class DataAnalysis(object):
    '''
 |   class DataAnalysis(object)
 |
 | This includes: 
 | (I)  XANES data processing
 |     1.1 Averaged/summed XANES plotting with interpolation for incident energy
 |     1.2 XANES area normalization (normalized to whole area or specified tail region)
 |     1.3 XANES area integration calculation
 |     1.4 Find peaks(maxima) in XANES
 |     1.5 Radiation damage XANES test
 | (II) RIXS data processing
 |     2.1 2D/3D RIXS plane plotting with interpolation for both incident energy and emission energy
 |         2.1.1 concentration correction
 |         2.1.2 IE versus EE plotting 
 |         2.1.3 IE versus ET plotting 
 |         2.1.4 RIXS plane CIE CEE CET cuts and integration plotting
 |     2.2 Averaging and merging for RIXS planes
 | (III)Save the data into RIXS txt files so that can be further used by other software, e.g, Matlab
 |
 |
 |  Methods
 |  ----------
 |
 |  -----------------------------------------
 |  ------------- XANES PART ----------------
 |  -----------------------------------------
 |
 |  XANES_data(): get XANES merged data ndarray from SPEC file
 |      return [incident energy, intensity], dtype = 1d ndarray
 |
 |  Radiation_damage(): Averaging for a step of XANES
 |      return [incident energy, intensity], dtype = 1d ndarray
 |
 |  XANES_normalize(): Normalize XANES to area into unity(whole area or specified tail area)
 |      return [incident energy, normalized_intensity], dtype = 1d ndarray
 |
 |  XANES_find_peaks(): Find XANES peaks and plotting
 |      return [peak energy, peak_intensity], dtype = 1d ndarray
 |
 |  XANES_area(): calculate XANES area for specified energy range
 |      return area, dtype = float
 |
 |
 |  -----------------------------------------
 |  ------------- RIXS PART -----------------
 |  -----------------------------------------
 |
 |  RIXS_data() : To get RIXS data ndarray from SPEC file
 |      return data ndarray [incident energy, emission energy, intensity]
 |
 |  RIXS_data_constantET() : To get RIXS data ndarray from SPEC file. In this type of scan, ET is fixed
 |                           This is different with RIXS_data where the emission energy is fixed.
 |      return data ndarray [incident energy, emission energy, intensity]
 |
 |  RIXS_merge(): To merge (sum up/average different RIXS data ndarray)
 |      return data ndarray [incident energy, emission energy, intensity]
 |
 |  RIXS_display() : To plot RIXS planes
 |      return None
 |
 |  RIXS_imshow() : To plot RIXS planes, using imshow instead of contour plotting to better present the raw data
 |      return None
 |
 |  RIXS_cut() : To do CIE, CET, CEE cuts
 |      return CIE, CET, CEE data ndarray [incident energy/energy transfer, intensity]
 |
 |  RIXS_integration() : Integration along incident energy and energy transfer
 |      return integrated data ndarray [incident energy/energy transfer, intensity]
 |
 |  -----------------------------------------
 |  ----------- GENERAL FUNCTIONS ----------- 
 |  ---------- This is not methods! --------- 
 |  -----------------------------------------
 |
 |  saveFile() : Save data into .dat file so that the data can be processed with other softwares
 |
 |  normalize_toArea() : Normalize XANES to area into unity(whole area or specified tail area)
 |      return [incident energy, normalized_intensity], dtype = 1d ndarray 
 |
 |  find_peaks(): Find XANES peaks and plotting
 |      return [peak energy, peak_intensity], dtype = 1d ndarray
 |
 |  XANES_area(): calculate XANES area for specified energy range
 |      return area, dtype = float
 |
 |
 |  Parameters
 |  ----------
 |  path : the filepath of Specfile
 |
    '''
    
    def __init__(self, path):
        self.path = path
        self.sf = SpecFile(path)
    
    def XANES_data(self, firstScan, lastScan, skipScan = [], interp_npt_1eV = 20, 
                   method = 'average', savetxt = False, channel = 'det_dtc'):
        """
        To get XANES merged data ndarray from SPEC file
        The incident energy for scans can be different

        Parameters
        ----------
        firstScan : the index of first scan, e.g, 71 corresponding to fscan '72.1'
        lastScan : the index of first scan
        skipScan : the problematic scans that you want to skip, 
                   e.g, [2,5,8] -- skip '3.1','6.1','9.1' scans 
        interp_npt_1eV : the number of interpolation points for 1 eV, default: 20 points for 1 eV
                         e.g, Incident Energy: 6535 eV - 6545 eV, 11 eV, 115 points, -----> 220 points
        method : 'average' or 'sum' for intensity
        channel : 'det_dtc' for HERFD-XAS, 'IF2' for conventional XAS
        savetxt: default True, save the ET, EE data as folders 
        Returns
        -------
        out : A 1d data ndarray [incident_Energy_interp, XANES_merge_inten]
              incident_Energy_interp -----> interpolated incident energy
              XANES_merge_inten -----> interpolated intensity
        """
        
        # Each scan has different incident energy points
        # this step finds the highest incident energy of the corresponding scans
        #             and the lowest incident energy
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
        
        # Define the scans list (skip the problematic scans)
        scanList = []
        for n in range(firstScan, lastScan + 1):
            if n not in skipScan:
                scanList.append(n)

        # We then fill the empty array with interpolated concentration corrected intensity
        # Use scipy.interpolate.interp1d to interpolate the points
        for n in range(firstScan, lastScan + 1):
            incident_Energy = self.sf[n].data_column_by_name('arr_hdh_ene')
            # corresponding intensity
            inten = self.sf[n].data_column_by_name(channel)/self.sf[n].data_column_by_name('I02') # Normalized to I02
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
    
    def Radiation_damage(self, firstScan, lastScan, scanStep, interp_npt_1eV = 20, 
                         method = 'average', savetxt = False, channel = 'det_dtc'):
        """
        To get XANES merged data ndarray for Radiation damage test
        The incident energy for scans can be different

        Parameters
        ----------
        firstScan : the index of first scan, e.g, 71 corresponding to fscan '72.1'
        lastScan : the index of first scan
        scanStep : the step for radiation damage
        interp_npt_1eV : the number of interpolation points for 1 eV, default: 20 points for 1 eV
                         e.g, Incident Energy: 6535 eV - 6545 eV, 11 eV, 115 points, -----> 220 points
        method : 'average' or 'sum' for intensity
        channel : 'det_dtc' for HERFD-XAS, 'IF2' for conventional XAS
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
        for n in range(firstScan, lastScan + 1, scanStep):
            incident_Energy = self.sf[n].data_column_by_name('arr_hdh_ene')
            # corresponding intensity
            inten = self.sf[n].data_column_by_name(channel)/self.sf[n].data_column_by_name('I02') # Normalized to I02
            # Define our interp1d function for incident energy, fill the outbound value with 0
            f_interp = interp.interp1d(incident_Energy, inten, bounds_error = False)
            # Find our interpolated incident energy and corresponding intensity
            incident_Energy_interp = np.linspace(incident_Energy_min, incident_Energy_max, incident_Energy_interp_npt)
            inten_interp = f_interp (incident_Energy_interp)
            # Fill the empty array with the interpolated data
            # Now we only interpolate in incident axis, not in emission axis
            XANES_inten_array[n-firstScan,:] = inten_interp
            print('adding the'+ str(n)+ ' scan')

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
        Normalize XANES to area into unity(whole area or specified tail area)
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
    
    def XANES_area(self, XANES_data, energy_range):
        """
        Calculate XANES area
        Parameters
        ----------
        XANES_data : The XANES_data output ndarray
        energy_range: (e1, e2)
                  e1: The starting energy for calculating the area, in eV
                  e2: The ending energy for calculating the area, in eV
        Returns
        -------
        out : XANES_area, dtype = float

        """
        # Calculate pre-edge area
        edge_area_index = np.where((XANES_data[0] >= energy_range[0]/1000) & (XANES_data[0] <= energy_range[1]/1000))
        edge_area_intensity = XANES_data[1][edge_area_index[0]]
        # Calculate the tail edge area
        edge_area = np.trapz(edge_area_intensity, dx=1)
        print('The edge area from %d eV to %d eV is :'%(energy_range[0], energy_range[1]) + str(edge_area) )
        return edge_area
    
    
    def XANES_find_peaks(self, XANES_data, energy_range = None, accuracy = (3,30), plot = True):
        """
        Find XANES peaks
        Parameters
        ----------
        XANES_data : the XANES_data or XANES_normalize output ndarray
        energy_range: default
        An energy number, e.g, Incident Energy: 6600 eV, 
                                                -----> will do normalization from 6600 eV to the end of tail feature
                                      Not difine-----> Normalization to the whole area
        Returns
        -------
        out : A data ndarray of the maxima [peak_energy, peak_intensity]

        """
        # Define the width of peaks that we want to detect
        peak_detect_width = np.arange(accuracy[0],accuracy[1])
        # Find peaks for the whole XANES_data range
        peak_index = signal.find_peaks_cwt(XANES_data[1],np.arange(accuracy[0],accuracy[1]))
        peak_dataList = np.array([XANES_data[0][peak_index],XANES_data[1][peak_index]])
        if energy_range == None:
            # plot the figures
            if plot == True:
                plt.plot(XANES_data[0]*1000,XANES_data[1])
                plt.scatter(peak_dataList[0]*1000,peak_dataList[1],c='r')
                plt.xlabel('Energy [eV]')
                plt.show()
            return peak_dataList
        else:
            e1 = energy_range[0]
            e2 = energy_range[1]
            # Define range interested
            # And pick out the specific peaks
            interested_index = np.where((peak_dataList[0]*1000 >= e1)&(peak_dataList[0]*1000 <= e2))
            range_peak_dataList = np.array([peak_dataList[0][interested_index],peak_dataList[1][interested_index]])

            # plot the figures
            if plot == True:
                plt.plot(XANES_data[0]*1000,XANES_data[1])
                plt.scatter(peak_dataList[0]*1000,peak_dataList[1],c='r')
                plt.xlim(e1,e2)
                y_max = max(range_peak_dataList[1])
                plt.ylim(0,y_max*1.2)
                plt.xlabel('Energy [eV]')
                plt.show()
            return range_peak_dataList

    def RIXS_data(self,firstScan, lastScan, concCorrecScan = False, interp_npt_1eV = 20, 
                  choice = 'EE', savetxt = False, unit = 'eV'):
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
        unit: Energy unit -----> 'eV' or 'KeV', default is 'eV' (in original Specfiles are in KeV)
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
        if concCorrecScan != False:
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
            if concCorrecScan == False:
                # don't do concentration correction for intensity
                correc_inten = self.sf[n].data_column_by_name('det_dtc')/(self.sf[n].data_column_by_name('I02'))
            else:
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
            if unit == 'eV':
                dataArray_EE[0] = dataArray_EE[0]*1000
                dataArray_EE[1] = dataArray_EE[1]*1000
                return dataArray_EE
            return dataArray_EE
        if choice == 'ET': 
            if unit == 'eV':
                dataArray_ET[0] = dataArray_ET[0]*1000
                dataArray_ET[1] = dataArray_ET[1]*1000
                return dataArray_ET
            return dataArray_ET
        
    def RIXS_data_constantET(self,firstScan, lastScan, concCorrecScan = False):

        incident_Energy = np.array([self.sf[i].data_column_by_name('mono.energy')[1] for i in range(firstScan, lastScan+1)]) 
        emission_Energy_firstScan = self.sf[firstScan].data_column_by_name('Spec.Energy')
        Energy_transfer = np.zeros_like(self.sf[firstScan].data_column_by_name('mono.energy'))

        # Fill the empty energy_transfer array
        for n in range(len(Energy_transfer)):
            Energy_transfer[n] = incident_Energy[0] - emission_Energy_firstScan[n]

        # Creat an empty array
        MDfci_correc_inten = np.zeros(((len(Energy_transfer), len(incident_Energy))))
        
        
        # We then fill the empty array with intensity
        for n in range(firstScan, lastScan + 1):
            # intensity for each scan at fixed incident energy
            
            if concCorrecScan == False:
                correc_inten = self.sf[n].data_column_by_name('det_dtc')/self.sf[n].data_column_by_name('I02') #Normalized to I02
            else:
                # Collect concentration correction intensity into an array
                concCorrec_inten = self.sf[concCorrecScan].data_column_by_name('det_dtc')/self.sf[concCorrecScan].data_column_by_name('I02') # Normalized to I02
                correc_inten = self.sf[n].data_column_by_name('det_dtc')/(self.sf[n].data_column_by_name('I02')*concCorrec_inten[n-firstScan]) 
    
            MDfci_correc_inten[:,n-firstScan] = correc_inten

        # Define Grids, EE_XX: incident energy array, EE_YY: emission energy array
        EE_XX, EE_YY = np.meshgrid(incident_Energy, Energy_transfer)
        RIXS_dataArray = np.array([EE_XX, EE_YY, MDfci_correc_inten])
        return RIXS_dataArray
        
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
    
    def RIXS_display(self, dataArray, title = 'RIXS',  choice = 'EE', mode = '2d',
                     savefig = False, normalize_to_Preedge = False,):
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
        normalize_to_Preedge : set pre-edge maximum as the max color(vmax), 
                               default = False: automatically choose the max of the whole peak max as vmax
        
        Returns
        -------
        out : ndarray
            Array of zeros with the given shape, dtype, and order.
    """
        if mode == '2d':
            # ------------ CONTOURF METHOD ------------
            # Transfer the energy from KeV scale into eV scale
            XX = dataArray[0]
            YY = dataArray[1]
            intensity = dataArray[2]
            plt.figure()
            plt.contourf(XX, YY, intensity, 20, cmap = cm.RdYlGn_r)
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
                if savefig == True:
                    plt.savefig(title+'.png',dpi=300)
                # plt.ylim(y_lim)
                # plt.xlim(x_lim)
                # x_lim = (6537,6544),y_lim=(5892,5902), 
            elif choice == 'ET':
                plt.axis('equal')
                plt.ylabel('Energy Transfer [eV]')
                if savefig == True:
                    plt.savefig(title+'.png',dpi=300)
        elif mode == '3d':
            from mpl_toolkits.mplot3d import Axes3D
            from matplotlib.ticker import LinearLocator, FormatStrFormatter
            XX = dataArray[0]
            YY = dataArray[1]
            intensity = dataArray[2]
            fig = plt.figure(figsize=(12,8))
            ax = fig.gca(projection='3d')
            plotting_3d = ax.plot_surface(XX, YY, intensity, 
                              cmap=cm.RdYlGn_r,
                              rstride=1, cstride=1,
                              linewidth=1, 
                              #antialiased=True,
                              vmin = 0,
                              vmax = 0.3,
                              alpha = 1
                             )
            # Customize the z axis.
            ax.zaxis.set_major_locator(LinearLocator(10))
            ax.zaxis.set_major_formatter(FormatStrFormatter('%.01f'))
            ax.xaxis.set_label('energy')
            ax.grid(False)
            # Add colorbar
            fig.colorbar(plotting_3d, shrink=0.5, aspect=5)
            

        return plt.show()
    
        
    def RIXS_imshow(self, dataArray):
        levels = np.linspace(dataArray[2].min(), dataArray[2].max(), 12)
        extent = (dataArray[0].min(), dataArray[0].max(), 
                  dataArray[1].min(), dataArray[1].max())

        plt.imshow(dataArray[2], extent=extent, 
                   origin='lower', aspect='auto',interpolation='nearest')
        plt.contour(dataArray[2], extent=extent, 
                    origin='lower', levels=levels,cmap=plt.cm.gray, linewidths=0.5)
        return plt.show()
    
    def RIXS_normalization(self, RIXS_data, XX_range =(), plot = False, 
                           levelnumber = 20, xlim = (6537.3,6544.8), ylim = (638.5, 647), 
                           savefig = False):
        """
        RIXS plane normalized to pre-edges (one needs to define pre-edge incident energy range)
        
        Parameters
        ----------
        RIXS_data : the RIXS_data [XX, YY, intensity]
        XX_range: default ()
                  A tuple of pre-edge range, e.g, XX_range =  (6538,6542)
                  -----> will do normalization to the maximum peak in the range of 6538 eV to 6542 eV
                  Default if the range not defined
                  -----> Normalization to the whole area
        Returns
        -------
        out : A data ndarray of the normalized RIXS data [RIXS_XX, RIXS_YY, RIXS_intensity]

        """

        # default incident_energy_range is the whole range
        if XX_range == ():
            XX_min = RIXS_data[0].min()
            XX_max = RIXS_data[0].max()
            XX_range=(XX_min,XX_max)

        # Crop the RIXS plane into pre-edge region by defining the incident_energy_range
        # This step is to ensure we are choosing the maximum of pre-edge without influence of main edge
        # First find indexes of this pre-edge region
        incident_E_index = np.where((RIXS_data[0][0,:] >= XX_range[0]) & 
                                    (RIXS_data[0][0,:] <= XX_range[1]))
        # Find the corresponding intensity of the pre-edge region
        crop_intensity = RIXS_data[2][:, incident_E_index[0]]

        # Find the maximum of pre-edge peak
        pre_edge_max = np.nanmax(crop_intensity)
        #print(pre_edge_max)

        # Normalization of RIXS_data intensity
        norm_intensity = RIXS_data[2]/pre_edge_max

        # Put new normalized data into a new data array
        norm_RIXS_dataArray = np.array([RIXS_data[0], 
                                        RIXS_data[1],
                                        norm_intensity])
        if plot == True:
            # Auto scale Plotting
            # Each contour have differenr intensity range, so we need different levels for contour plotting
            intensity_range = np.nanmax(norm_intensity) - np.nanmin(norm_intensity)
            #level = int(intensity_range * levelscale)
            level = list(np.arange(0,1,1/levelnumber))
            # print('intensity range', intensity_range)
            fig = plt.figure(figsize=(6,6))
            MyContour = plt.contourf(norm_RIXS_dataArray[0],
                                     norm_RIXS_dataArray[1],
                                     norm_RIXS_dataArray[2],
                                     levels = level,
                                     cmap=cm.RdYlGn_r,
                                     vmin = -0.2,
                                     vmax = 1,
                                     extend="both",
                                    )
            MyContour.cmap.set_over('#aa1b24')
            MyContour.cmap.set_under('white')
            plt.colorbar()
            plt.contour(norm_RIXS_dataArray[0],
                        norm_RIXS_dataArray[1],
                        norm_RIXS_dataArray[2],
                        levels = level,
                        linewidths = 0.3, 
                        colors='black')
            plt.xlim(xlim[0], xlim[1])
            plt.ylim(ylim[0], ylim[1])
            #plt.yticks(np.arange(640, 648, 1))
            plt.gca().set_aspect('equal', adjustable='box')
            # plt.title(title)
            plt.xlabel('Incident Energy [eV]')
            plt.ylabel('Energy Transfer [eV]')
            if savefig == True:
                fig.savefig('norm_RIXS', dpi =300)
            plt.show()
            
        return norm_RIXS_dataArray
    
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
    
    def RIXS_integration(self, dataArray, choice = 'IE'):
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
        if choice == 'IE':
            integration_dataArray = np.array([dataArray[0][0,:]*1000,sumIntensity_IE])
        elif choice == 'ET':
            integration_dataArray = np.array([dataArray[1][:,0]*1000,sumIntensity_ET])
        #integration_dataArray = np.array([[dataArray[0][0,:]*1000,sumIntensity_IE],[dataArray[1][:,0]*1000,sumIntensity_ET]])
        return integration_dataArray

    
# Functions
def saveFile(dataList, headerList, folderPath = 'None', fileName = 'myData', choice = 'XANES'):
    """
    Save data into txt
    Parameters
    ----------
    dataList : [dataArray1, dataArray2..., dataArrayn]
    headerList : ['Column title 1', 'Column title 2'..., 'Column title n']
    folderPath: the folder where you want to save your data 
    fileName: the file name
    choice: 'XANES' or 'RIXS'
    
    Returns
    -------
    out : A '.dat' file saving all the data

    """
    if folderPath == 'None':
        folderPath = os.getcwd()
    
    headerSaveList = ' , '.join(headerList)
    # Save XANES
    if choice == 'XANES':
        dataSaveList = np.array(dataList) 
        np.savetxt(folderPath + fileName + '.dat', 
                   np.transpose(dataSaveList), fmt = '%.12f', 
                   header = headerSaveList)
    # Save RIXS
    elif choice == 'RIXS':
        np.savetxt(folderPath + fileName + '.dat', 
                   dataList, fmt = '%.12f', 
                   header = headerSaveList)
    
def normalize_toArea(XANES_data, normalized_starting_energy = None):
    """
    Normalize XANES to area into unity(whole area or specified tail area)
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
    norm_intensity = 20*XANES_data[1]/tail_edge_area
    norm_dataArray = np.array([XANES_data[0],norm_intensity])
    return norm_dataArray
 
def find_peaks(XANES_data, energy_range = None, accuracy = (3,30), plot = True):
    """
    Find XANES peaks
    Parameters
    ----------
    XANES_data : the XANES_data or XANES_normalize output ndarray
    energy_range: default
    An energy number, e.g, Incident Energy: 6600 eV, 
                                          -----> will do normalization from 6600 eV to the end of tail feature
                                Not difine-----> Normalization to the whole area
    Returns
    -------
    out : A data ndarray of the maxima [peak_energy, peak_intensity]

    """
    # Define the width of peaks that we want to detect
    peak_detect_width = np.arange(accuracy[0],accuracy[1])
    # Find peaks for the whole XANES_data range
    peak_index = signal.find_peaks_cwt(XANES_data[1],np.arange(accuracy[0],accuracy[1]))
    peak_dataList = np.array([XANES_data[0][peak_index],XANES_data[1][peak_index]])
    if energy_range == None:
        # plot the figures
        if plot == True:
            plt.plot(XANES_data[0]*1000,XANES_data[1])
            plt.scatter(peak_dataList[0]*1000,peak_dataList[1],c='r')
            plt.xlabel('Energy [eV]')
            plt.show()
        return peak_dataList
    else:
        e1 = energy_range[0]
        e2 = energy_range[1]
        # Define range interested
        # And pick out the specific peaks
        interested_index = np.where((peak_dataList[0]*1000 >= e1)&(peak_dataList[0]*1000 <= e2))
        range_peak_dataList = np.array([peak_dataList[0][interested_index],peak_dataList[1][interested_index]])

        # plot the figures
        if plot == True:
            plt.plot(XANES_data[0]*1000,XANES_data[1])
            plt.scatter(peak_dataList[0]*1000,peak_dataList[1],c='r')
            plt.xlim(e1,e2)
            y_max = max(range_peak_dataList[1])
            plt.ylim(0,y_max*1.2)
            plt.xlabel('Energy [eV]')
            plt.show()
        return range_peak_dataList


def XANES_area(XANES_data, energy_range):
    """
    Calculate XANES area
    Parameters
    ----------
    XANES_data : The XANES_data output ndarray
    energy_range: (e1, e2)
              e1: The starting energy for calculating the area, in eV
              e2: The ending energy for calculating the area, in eV
    Returns
    -------
    out : XANES_area, dtype = float

    """
    # Calculate pre-edge area
    edge_area_index = np.where((XANES_data[0] >= energy_range[0]/1000) & (XANES_data[0] <= energy_range[1]/1000))
    edge_area_intensity = XANES_data[1][edge_area_index[0]]
    # Calculate the tail edge area
    edge_area = np.trapz(edge_area_intensity, dx=1)
    print('The edge area from %d eV to %d eV is :'%(energy_range[0], energy_range[1]) + str(edge_area) )
    return edge_area
