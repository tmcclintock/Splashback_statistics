import numpy as np
import pandas as pd
from colossus.cosmology import cosmology
from colossus.lss import peaks

#Header of the catalogs. don't touch this :)
_default_header = "#scale(0) id(1) desc_scale(2) desc_id(3) num_prog(4) pid(5) upid(6) desc_pid(7) phantom(8) sam_M200b(9) M200b(10) R200b(11) rs(12) vrms(13) mmp?(14) scale_of_last_MM(15) vmax(16) x(17) y(18) z(19) vx(20) vy(21) vz(22) Jx(23) Jy(24) Jz(25) Spin(26) Breadth_first_ID(27) Depth_first_ID(28) Tree_root_ID(29) Orig_halo_ID(30) Snap_num(31) Next_coprogenitor_depthfirst_ID(32) Last_progenitor_depthfirst_ID(33) Last_mainleaf_depthfirst_ID(34) Tidal_Force(35) Tidal_ID(36) Rs_Klypin(37) M200b_all(38) Mvir(39) M200c(40) M500c(41) M2500c(42) Xoff(43) Voff(44) Spin_Bullock(45) b_to_a(46) c_to_a(47) A[x](48) A[y](49) A[z](50) b_to_a_500c(500c)(51) c_to_a_500c(500c)(52) A[x]_500c(500c)(53) A[y]_500c(500c)(54) A[z]_500c(500c)(55) T/|U|(56) M_pe_Behroozi(57) M_pe_Diemer(58) Halfmass_Radius(59) Macc(60) Mpeak(61) Vacc(62) Vpeak(63) Halfmass_Scale(64) Acc_Rate_Inst(65) Acc_Rate_100Myr(66) Acc_Rate_1*Tdyn(67) Acc_Rate_2*Tdyn(68) Acc_Rate_Mpeak(69) Mpeak_Scale(70) Acc_Scale(71) First_Acc_Scale(72) First_Acc_Mvir(73) First_Acc_Vmax(74) Vmax\@Mpeak(75) Tidal_Force_Tdyn(76) Log_[Vmax/Vmax_max[Tdyn;Tmpeak]](77) Time_to_future_merger(78) Future_merger_MMP_ID(79) Rsp_status(80) upid_mean(81) Rsp_mean(82) Rsp_mean_err(83) Msp_mean(84) Msp_mean_err(85) upid_percentile50(86) Rsp_percentile50(87) Rsp_percentile50_err(88) Msp_percentile50(89) Msp_percentile50_err(90) upid_percentile75(91) Rsp_percentile75(92) Rsp_percentile75_err(93) Msp_percentile75(94) Msp_percentile75_err(95) upid_percentile87(96) Rsp_percentile87(97) Rsp_percentile87_err(98) Msp_percentile87(99) Msp_percentile87_err(100)"

class Catalog(object):
    """Splashback halo catalog
    
    Args:
        length (float): side length of the simulation
        scale_factor (float): scale factor of the snapshot
        cosmo (string): cosmology of the sim. Default is bolshoi
    """
    def __init__(self, length, scale_factor, cosmo="bolshoi"):
        names    = _default_header.split(" ")
        for i,name in enumerate(names):
            newname = name[:name.find('(')]
            names[i] = newname
        names[0] = names[0][1:] #remove the pound sign from scale

        #Load in the catalog as an array
        data = np.loadtxt("L%04d_N1024_CBol/hlist_%.5f_mpeak.list"%(length, scale_factor))

        #Create the dataframe
        df = pd.DataFrame(data=data, columns=names)
        
        #Add columns for the things that we want
        #This includes c200b and the ratios X_
        df["c200b"] = df["R200b"].values/df["rs"].values
        for kind in ["mean", "mean_err", "percentile50",
                     "percentile75", "percentile87"]:
            df["X_Rsp_%s"%kind] = df["Rsp_%s"%kind].values / df["R200b"]
            df["X_Msp_%s"%kind] = df["Msp_%s"%kind].values / df["M200b"]
            continue
        #Now the peak heights
        cosmology.setCosmology(cosmo)
        redshift = 1./scale_factor - 1
        for kind in ["200b", "sp_mean", "sp_percentile50",
                     "sp_percentile75", "sp_percentile87"]:
            df["nu%s"%kind] = peaks.peakHeight(df["M%s"%kind].values, redshift)
            continue
        #Make the dataframe an attribute and that's it
        self.dataframe = df

        #Make a dictionary of the correlated variables
        self.correlated_variables = {}

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

    def compute_correlations(self, R_or_M, kind="mean"):
        """Compute the correlations between either X_R or X_M
        with everything else.

        Stores the correlation attributes in a dictionary, and
        also returns them.

        Args:
            R_or_M (string): either 'R' or 'M'
            kind (string): the kind of splashback definition, e.g. "mean" or "percentile75"

        Returns:
            names: array of strings ordered by absolute correlation.
            correlations: array of correlations

        """
        if R_or_M not in ["R", "M"]:
            raise Exception("R_or_M must be either 'R' or 'M'.")
        if kind not in ["mean", "percentile50",
                        "percentile75", "percentile87"]:
            raise Exception("kind must be 'mean', 'percentile50', 'percentile75', or 'percentile87'.")
        
        if "%s_%s"%(R_or_M, kind) in self.correlated_variables.keys():
            #Correlations already computed
            return

        X = self.dataframe['X_%ssp_%s'%(R_or_M, kind)].values
        ordered_names = np.array([])
        ordered_corrs = np.array([])
        for name in self.dataframe.columns:
            v = self.dataframe[name].values
            ordered_corrs = np.append(ordered_corrs, np.corrcoef(X, v)[0,1])
            ordered_names = np.append(ordered_names, name)
        order = np.argsort(np.fabs(ordered_corrs))
        self.ordered_corrs = ordered_corrs[order[::-1]]
        self.ordered_names = ordered_names[order[::-1]]

        #Save the correlations
        self.correlated_variables['%s_%s'%(R_or_M, kind)] = \
            [self.ordered_names, self.ordered_corrs]

        #Return the result
        return self.ordered_names, self.ordered_corrs
    

if __name__ == "__main__":
    cat = Catalog(2000, 1.)
    print(cat.dataframe.columns)
    names, corrs = cat.compute_correlations("R")

    for name, R in zip(names, corrs):
        print(name, R)
