import numpy as np
import pandas as pd
from colossus.cosmology import cosmology
from colossus.lss import peaks
from minepy import MINE

#Header of the catalogs. don't touch this :)
_default_header = "#scale(0) id(1) desc_scale(2) desc_id(3) num_prog(4) pid(5) upid(6) desc_pid(7) phantom(8) sam_M200b(9) M200b(10) R200b(11) rs(12) vrms(13) mmp?(14) scale_of_last_MM(15) vmax(16) x(17) y(18) z(19) vx(20) vy(21) vz(22) Jx(23) Jy(24) Jz(25) Spin(26) Breadth_first_ID(27) Depth_first_ID(28) Tree_root_ID(29) Orig_halo_ID(30) Snap_num(31) Next_coprogenitor_depthfirst_ID(32) Last_progenitor_depthfirst_ID(33) Last_mainleaf_depthfirst_ID(34) Tidal_Force(35) Tidal_ID(36) Rs_Klypin(37) M200b_all(38) Mvir(39) M200c(40) M500c(41) M2500c(42) Xoff(43) Voff(44) Spin_Bullock(45) b_to_a(46) c_to_a(47) A[x](48) A[y](49) A[z](50) b_to_a_500c(500c)(51) c_to_a_500c(500c)(52) A[x]_500c(500c)(53) A[y]_500c(500c)(54) A[z]_500c(500c)(55) T/|U|(56) M_pe_Behroozi(57) M_pe_Diemer(58) Halfmass_Radius(59) Macc(60) Mpeak(61) Vacc(62) Vpeak(63) Halfmass_Scale(64) Acc_Rate_Inst(65) Acc_Rate_100Myr(66) Acc_Rate_1*Tdyn(67) Acc_Rate_2*Tdyn(68) Acc_Rate_Mpeak(69) Mpeak_Scale(70) Acc_Scale(71) First_Acc_Scale(72) First_Acc_Mvir(73) First_Acc_Vmax(74) Vmax\@Mpeak(75) Tidal_Force_Tdyn(76) Log_[Vmax/Vmax_max[Tdyn;Tmpeak]](77) Time_to_future_merger(78) Future_merger_MMP_ID(79) Rsp_status(80) upid_mean(81) Rsp_mean(82) Rsp_mean_err(83) Msp_mean(84) Msp_mean_err(85) upid_percentile50(86) Rsp_percentile50(87) Rsp_percentile50_err(88) Msp_percentile50(89) Msp_percentile50_err(90) upid_percentile75(91) Rsp_percentile75(92) Rsp_percentile75_err(93) Msp_percentile75(94) Msp_percentile75_err(95) upid_percentile87(96) Rsp_percentile87(97) Rsp_percentile87_err(98) Msp_percentile87(99) Msp_percentile87_err(100)"

skipped_names = ["scale","id","desc_scale","desc_id","num_prog","pid","upid",
                 "phantom","sam_M200b","rs","vx","vy","vz","mmp?","x","y","z","Jx",
                 "Jy","Jz","Breadth_first_ID","Depth_first_ID",
                 "Last_mainleaf_depthfirst_ID","Tidal_ID","Tree_root_ID",
                 "Orig_halo_ID","Snap_num","Next_coprogenitor_depthfirst_ID",
                 "Last_progenitor_depthfirst_ID","Last_mainleaf_depthfirst_ID",
		 "c_to_a","b_to_a","c_to_a_500c","b_to_a_500c",
                 "M200b","R200b","M200b_all","Mvir","M200c","M500c","M2500c",
                 "M_pe_Behroozi","M_pe_Diemer","Macc","Mpeak","First_Acc_Mvir",
                 "First_Acc_Vmax","Vmax\@Mpeak","Rsp_status",
                 "Future_merger_MMP_ID","Time_to_future_merger","desc_pid"]

class Catalog(object):
    """Splashback halo catalog
    
    Args:
        lengths (array-like): list of side lengths of the simulations
        scale_factor (float): scale factor of the snapshot
        cosmo (string): cosmology of the sim. Default is bolshoi
    """
    def __init__(self, lengths, scale_factor, cosmo="bolshoi", parents_only=True):
        names    = _default_header.split(" ")
        for i, name in enumerate(names):
            newname = name[:name.find('(')]
            names[i] = newname
        names[0] = names[0][1:] #remove the pound sign from scale

	# if the number of input catalog = 1
        if len(lengths) == 1:
            #Load in the catalog as an array
            data = np.loadtxt("sparta_cats/L%04d_N1024_CBol/hlist_%.5f_mpeak.list"%(lengths[0], scale_factor))

            #Create the dataframe
            df = pd.DataFrame(data=data, columns=names)
        #if we want to concatenate all the input catalogs
        else:
            for i in range(len(lengths)):
                if i == 0:
                    df = pd.DataFrame(data = np.loadtxt("sparta_cats/L%04d_N1024_CBol/hlist_%.5f_mpeak.list"%(lengths[i], scale_factor)), columns=names)
                else:
                    df = pd.concat([df,pd.DataFrame(data = np.loadtxt("sparta_cats/L%04d_N1024_CBol/hlist_%.5f_mpeak.list"%(lengths[i], scale_factor)), columns=names)])

        #Add columns for the things that we want
        #This includes c200b and the ratios X_
        df["c200b"] = df["R200b"].values/df["rs"].values
        for kind in ["mean", "mean_err", "percentile50",
                     "percentile75", "percentile87"]:
            df["X_Rsp_%s"%kind] = df["Rsp_%s"%kind].values / df["R200b"]
            df["X_Msp_%s"%kind] = df["Msp_%s"%kind].values / df["M200b"]
            continue
        # Define 3d ellipticity
        df["e_3d"] = (df["c_to_a"].values - 1.)/(2.*(df["c_to_a"].values+df["b_to_a"].values+1.))
        df["e_3d_500c"] = (df["c_to_a_500c"].values - 1.)/(2.*(df["c_to_a_500c"].values+df["b_to_a_500c"].values+1.))

        # Normalize the radius-related property.
        df["Halfmass_Radius"] /= df["R200b"]

        # Normalize the accretion rate
        # for now, we assume we are working on a single snapshot, but if we want to generalize it to multiple snapshots, we have to remove the scale-factor dependency too.
        # Acc_rockstar = dM/dt = (dM/dlnM)*(dlnM/dlna)*(dlna/da)*(da/dt)
        # Thus, dlnM/dlna = (Acc_rockstar/M)*(a*dt/da) -- note that (a*dt/da) is consistent within the same snapshot		
        df["Acc_Rate_Inst"] /= df["Msp_mean"]
        df["Acc_Rate_100Myr"] /= df["Msp_mean"]
        df["Acc_Rate_1*Tdyn"] /= df["Msp_mean"]
        df["Acc_Rate_2*Tdyn"] /= df["Msp_mean"]
        df["Acc_Rate_Mpeak"] /= df["Msp_mean"]

        #Now the peak heights
        cosmology.setCosmology(cosmo)
        redshift = 1./scale_factor - 1
        for kind in ["200b", "sp_mean", "sp_percentile50",
                     "sp_percentile75", "sp_percentile87"]:
            df["nu%s"%kind] = peaks.peakHeight(df["M%s"%kind].values, redshift)
            continue

        #Make the dataframe an attribute and that's it
        self.dataframe = df
        #Make a dictionary of the correlated variables and MICed variables
        self.correlated_variables = {}
        self.MICed_variables = {}

    def property(self, name):
        """Get the values of a property in the dataframe.

        Args:
            name (string): column name

        Return:
            column: numpy array of the values in the halo catalog
        """
        if name not in self.dataframe.columns:
            raise Exception("Invalid property.")
        return self.dataframe[name].values

    def add_property(self, name, values):
        self.dataframe[name] = values
        return

    def compute_correlations(self, R_or_M, kind="mean"):
        """Compute the correlations between either X_R or X_M
        with everything else.

        Stores the correlation attributes in a dictionary, and
        also returns them.

        Args:
            R_or_M (string): either 'R' or 'M'
            kind (string): the kind of splashback definition, e.g. 
                "mean" or "percentile75"

        Returns:
            names: array of strings ordered by absolute correlation.
            correlations: array of correlations

        """
        if R_or_M not in ["R", "M"]:
            raise Exception("R_or_M must be either 'R' or 'M'.")
        if kind not in ["mean", "percentile50",
                        "percentile75", "percentile87"]:
            raise Exception("kind must be 'mean', 'percentile50', "+
                            "'percentile75', or 'percentile87'.")

        this_var = "%s_%s"%(R_or_M, kind)
        if this_var in self.correlated_variables.keys():
            #Correlations already computed
            return

        pids = self.dataframe['upid_%s'%kind] #pics out parent halos
        inds = pids < 0
        
        X = self.dataframe['X_%ssp_%s'%(R_or_M, kind)].values
        ordered_names = np.array([])
        ordered_corrs = np.array([])
        for name in self.dataframe.columns:
            #Skip certain keywords
            if name in skipped_names:
                continue
            if any(x in name for x in ["upid","percentile",
                                       'X_%ssp_%s'%(R_or_M, kind)]):
                continue
            v = self.dataframe[name].values
            ordered_corrs = np.append(ordered_corrs, np.corrcoef(X[inds], v[inds])[0,1])
            ordered_names = np.append(ordered_names, name)
        order = np.argsort(np.fabs(ordered_corrs))
        #Ascending order
        self.ordered_corrs = ordered_corrs[order][::-1]
        self.ordered_names = ordered_names[order][::-1]
        self.order = order

        #Save the correlations
        self.correlated_variables['%s_%s'%(R_or_M, kind)] = \
            [self.ordered_names, self.ordered_corrs]

        #Return the result
        return self.ordered_names, self.ordered_corrs

    def compute_partial_correlations(self, R_or_M, y, kind="mean"):
        """Compute the correlations between either X_R or X_M
        with everything else, with any one variable (y) FIXED.

        Stores the partial correlation attributes in a dictionary, and
        also returns them.

        Args:
            R_or_M (string): either 'R' or 'M'
            y (string): the name of the "fixed" variable
            kind (string): the kind of splashback definition, e.g. 
                "mean" or "percentile75"

        Returns:
            names: array of strings ordered by absolute partial correlation.
            correlations: array of partial correlations

        """
        if R_or_M not in ["R", "M"]:
            raise Exception("R_or_M must be either 'R' or 'M'.")
        if kind not in ["mean", "percentile50",
                        "percentile75", "percentile87"]:
            raise Exception("kind must be 'mean', 'percentile50', "+
                            "'percentile75', or 'percentile87'.")

        pids = self.dataframe['upid_%s'%kind] #pics out parent halos
        inds = pids < 0

        y_vec = self.dataframe[y].values  ### the fixed variable
        X = self.dataframe['X_%ssp_%s'%(R_or_M, kind)].values
        names, two_point_corr = self.correlated_variables['%s_%s'%(R_or_M, kind)] ### correlations should have been calculated before running this function
        corr_Xy = np.corrcoef(X[inds],y_vec[inds])[0,1]  #dependent variable (X) vs fixed variable
        ordered_names = np.array([])
        ordered_corrs = np.array([])
        for name in self.dataframe.columns:
            #Skip certain keywords
            if name in skipped_names:
                continue
            if name == y:
                continue
            if any(x in name for x in ["upid","percentile",
                                       'X_%ssp_%s'%(R_or_M, kind)]):
                continue
            idx_var = np.where(names==name)[0]
            corr_Xv = two_point_corr[idx_var] #dependent variable (X) vs variable of interest
            v = self.dataframe[name].values[inds]
            corr_vy = np.corrcoef(v,y_vec[inds])[0,1] #fixed variable vs variable of interest
            #now calculate partial correlation 
            ordered_corrs = np.append(ordered_corrs, (corr_Xv-corr_Xy*corr_vy)/np.sqrt(1.-corr_Xy**2)*np.sqrt(1.-corr_vy**2))
            ordered_names = np.append(ordered_names, name)	

        order = np.argsort(np.fabs(ordered_corrs))
        ordered_corrs_partial = ordered_corrs[order][::-1]
        ordered_names_partial = ordered_names[order][::-1]

        self.correlated_variables['%s_%s_%s-fixed'%(R_or_M, kind,y)] = \
            [ordered_names_partial, ordered_corrs_partial]

        return ordered_names_partial, ordered_corrs_partial

    def compute_multiple_correlations(self, R_or_M, variables, kind="mean"):
        """Compute the multiple correlation between either X_R or X_M
        with A SET OF variables.

        Args:
            R_or_M (string): either 'R' or 'M'
            variables (array of string): the independent variables of interest
            kind (string): the kind of splashback definition, e.g. 
                "mean" or "percentile75"

        Returns:
            names: array of strings ordered by absolute partial correlation.
            correlations: array of partial correlations

        """
        if R_or_M not in ["R", "M"]:
            raise Exception("R_or_M must be either 'R' or 'M'.")
        if kind not in ["mean", "percentile50",
                        "percentile75", "percentile87"]:
            raise Exception("kind must be 'mean', 'percentile50', "+
                            "'percentile75', or 'percentile87'.")
        pids = self.dataframe['upid_%s'%kind] #pics out parent halos
        inds = pids < 0

        names, two_point_corr = self.correlated_variables['%s_%s'%(R_or_M, kind)]
        ind_name = np.in1d(names,variables)
        name_ordered = names[ind_name]
        corr_Xv = two_point_corr[ind_name]

        corr_matrix_of_vars = np.zeros((len(variables),len(variables)))
        for i in range(len(variables)):
            for j in range(len(variables)):
                arr1 = (self.dataframe[name_ordered[i]].values)[inds]
                arr2 = (self.dataframe[name_ordered[j]].values)[inds]
                corr_matrix_of_vars[i,j] = np.corrcoef(arr1,arr2)[0,1]

        multiple_corr = np.matmul(np.matmul(corr_Xv,np.linalg.inv(corr_matrix_of_vars)),corr_Xv.T)
        return name_ordered, np.sqrt(multiple_corr)

    def compute_MICs(self, R_or_M, kind="mean",
                     alpha=0.6, c=15, est="mic_approx"):
        """Compute the maximal information coefficient (MIC) between either X_R or X_M
        with everything else.

        Stores the MIC attributes in a dictionary, and
        also returns them.

        Args:
            R_or_M (string): either 'R' or 'M'
            kind (string): the kind of splashback definition, e.g. 
                "mean" or "percentile75"
            alpha (float): see MINEPY docs
            c (int): see MINEPY docs
            est (string): see MINEPY docs

        Returns:
            names: array of strings ordered by MIC values.
            MICs: array of MICs

        """
        if R_or_M not in ["R", "M"]:
            raise Exception("R_or_M must be either 'R' or 'M'.")
        if kind not in ["mean", "percentile50",
                        "percentile75", "percentile87"]:
            raise Exception("kind must be 'mean', 'percentile50', "+
                            "'percentile75', or 'percentile87'.")
        
        if "%s_%s"%(R_or_M, kind) in self.MICed_variables.keys():
            #MICs already computed
            return
        
        pids = self.dataframe['upid_%s'%kind] #pics out parent halos
        inds = pids < 0

        #Make a MINE object to computes MICs
        mine = MINE(alpha=0.6, c=15, est="mic_approx")
        
        X = self.dataframe['X_%ssp_%s'%(R_or_M, kind)].values
        MIC_ordered_names = np.array([])
        ordered_MICs = np.array([])
        for name in self.dataframe.columns:
            if name in skipped_names:
                continue
            if any(x in name for x in ["upid","percentile",
                                       'X_%ssp_%s'%(R_or_M, kind)]):
                continue
            v = self.dataframe[name].values
            mine.compute_score(X[inds], v[inds])
            ordered_MICs = np.append(ordered_MICs, mine.mic())
            MIC_ordered_names = np.append(MIC_ordered_names, name)
        MIC_order = np.argsort(ordered_MICs)
        #Ascending order
        self.ordered_MICs = ordered_MICs[MIC_order][::-1]
        self.MIC_ordered_names = MIC_ordered_names[MIC_order][::-1]
        self.MIC_order = MIC_order

        #Save the MICs
        self.MICed_variables['%s_%s'%(R_or_M, kind)] = \
            [self.MIC_ordered_names, self.ordered_MICs]

        #Return the result
        return self.MIC_ordered_names, self.ordered_MICs

        
if __name__ == "__main__":
#    cat = Catalog([2000,1000,500,250,125,63], 0.88534)
    cat = Catalog([2000], 0.88534)
    names, corrs = cat.compute_correlations("R", "mean")
    names_partial_acc, corrs_partial_acc = cat.compute_partial_correlations("R","Acc_Rate_2*Tdyn","mean")
    names_partial_nu, corrs_partial_nu = cat.compute_partial_correlations("R","nusp_mean","mean")
    print("-----------------------------------")

    print("###################################")
    print("regular correlation between X_Rsp_mean and variables")
    print("###################################")
    for i in range(len(names)):
        print("%20s %.4f"%(names[i], corrs[i]))
    print("-----------------------------------")

    print("###################################")
    print("partial correlation between X_Rsp_mean and variables, FIXING Acc_Rate_2*Tdyn")
    print("###################################")
    for i in range(len(names_partial_acc)):
        print("%20s %.4f"%(names_partial_acc[i], corrs_partial_acc[i]))
    print("-----------------------------------")

    print("###################################")
    print("partial correlation between X_Rsp_mean and variables, FIXING nusp_mean")
    print("###################################")
    for i in range(len(names_partial_nu)):
        print("%20s %.4f"%(names_partial_nu[i], corrs_partial_nu[i]))
    
    print("-----------------------------------")
    print(cat.compute_multiple_correlations("R",["Acc_Rate_2*Tdyn","nusp_mean"],"mean"))

    print(cat.compute_multiple_correlations("R",["Acc_Rate_2*Tdyn","T/|U|"],"mean"))
    print(cat.compute_multiple_correlations("R",["T/|U|","nusp_mean"],"mean"))

    print(cat.compute_multiple_correlations("R",["Acc_Rate_2*Tdyn","Halfmass_Scale"],"mean"))
    print(cat.compute_multiple_correlations("R",["nusp_mean","Halfmass_Scale"],"mean"))

    print(cat.compute_multiple_correlations("R",["Acc_Rate_2*Tdyn","Spin"],"mean"))
    print(cat.compute_multiple_correlations("R",["nusp_mean","Spin"],"mean"))
