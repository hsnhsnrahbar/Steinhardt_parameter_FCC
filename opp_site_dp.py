import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import math
import pandas as pd
import openpyxl
import cantera as ct
import os
#import ch_sol



# import chemical_kinetics_solver  # Commented out as it's not a standard library

# Uncomment this line if you have the chemical_kinetics_solver module
# gas = ct.Solution('FFCM2.yaml') 
##
# Constants and parameters
H_diameter = 66e-12
C_diameter = 154e-12
O_diameter = 146e-12
N_diameter = 148e-12

molecule_area_c2h2 = (3.14*2*(H_diameter**2) + 3.14*2*(C_diameter**2))/4
molecule_area_c4h4 = (3.14*4*(H_diameter**2) + 3.14*4*(C_diameter**2))/4

# Constants
NAVA = 6.022e23  # Avogadro's number
kb = 1.38065e-23  # Boltzmann constant, J/K

# 
gas = ct.Solution('FFCM2.yaml') 
kinetics_mechanism_type="FFCM2" 

xdtube = .054 # Burner tube diameter (m)
pressure_atm=101325/101325

vel = .001
reactor_length=.2 #in m


# Ma et al. method for growth rate calculation
class MaGrowthModel:
    def __init__(self, dp=15e-9, P_FEEDSTOCK_c2h2=20000.645, P_FEEDSTOCK_c4h4=20000.645):
        self.dp = dp
        self.P_FEEDSTOCK_c2h2 = P_FEEDSTOCK_c2h2
        self.P_FEEDSTOCK_c4h4 = P_FEEDSTOCK_c4h4
        self.NAVA = 6.022e23
        self.kb = 1.38065e-23
        self.R_gas = 8.314
        
        # Atomic carbon volume (m³/atom)
        self.volume_carbon_atom = 1.44e-30
        
        # Gas properties
        self.mass_C2H2 = 4.324e-26  # kg/molecule
        self.mass_C4H4 = 8.66e-26 # kg/molecule
        
        # Constants from Ma et al.
        self.D0 = 0.5e-4  # m²/s 
        self.E_act = 1.5 * 1.602e-19  # J (1.5 eV)
        
    def calculate_gamma(self, temp_K):
        """Calculate reaction probability based on temperature"""
        temp_C = temp_K - 273.15
        
        if temp_C <= 700:
            gamma = 6.382123e-03 * np.exp(1.113002e-02 * (temp_K - 823.81))
        else:
            gamma = 3.461503e-02 * np.exp(-3.268908e-02 * (temp_K - 973.00))
        return gamma
    
    def calculate_growth_rate_per_step(self, temp_K):
        """Calculate CNT growth rate in m/s using Ma et al. method"""
        # Site density
        particle_area = np.pi * (self.dp)**2/2# /2 is due to the fact that half of the NP is inactive since cnt covres that
        
        '''
        which calculation should be used to measure the number of sites
        '''
       # site_per_particle = particle_area/molecule_area_c2h2/2/1.8
        site_per_particle = particle_area/((self.P_FEEDSTOCK_c2h2/(self.P_FEEDSTOCK_c2h2+self.P_FEEDSTOCK_c4h4))*molecule_area_c2h2+(self.P_FEEDSTOCK_c4h4/(self.P_FEEDSTOCK_c2h2+self.P_FEEDSTOCK_c4h4))*molecule_area_c4h4)/1.8/2




        rho_particle = 1e15
        rho_0 = site_per_particle * rho_particle
        
        # Gas concentration
        n_c2h2 = self.NAVA * self.P_FEEDSTOCK_c2h2 / (self.R_gas * temp_K)
        n_c4h4 = self.NAVA * self.P_FEEDSTOCK_c4h4 / (self.R_gas * temp_K)
        
        # Reaction probability
        gamma = self.calculate_gamma(temp_K)
        projected_area_per_particle = np.pi * (self.dp/2)**2
        
        # Surface coverage factor
        A_star = rho_particle * projected_area_per_particle
   
        # Forward rate constants
        v_bar1 = math.sqrt(self.R_gas * temp_K / (2 * math.pi * self.mass_C2H2 * self.NAVA))
        kf1 = v_bar1 * gamma * (A_star / rho_0)
        v_bar2 = math.sqrt(self.R_gas * temp_K / (2 * math.pi * self.mass_C4H4 * self.NAVA))
        kf2 = v_bar2 * gamma * (A_star / rho_0)
        
        # Precipitation rate constant
        R = self.dp/2
        k3 = (self.D0 / R**2) * math.exp(-self.E_act / (self.kb * temp_K))
        
        # Steady-state active sites
        O_ss = (k3 * rho_0) / ((kf1 * n_c2h2 + kf2 * n_c4h4) + k3)
        
        # Reaction rates
        q_C2H2 = kf1 * n_c2h2 * O_ss
        q_C4H4 = kf2 * n_c4h4 * O_ss
        
        # Carbon atoms flux
        carbon_atoms_flux = (2 * q_C2H2 + 4 * q_C4H4)
        
        # Growth rate calculation
        d_outer = self.dp
        d_inner = d_outer/2
        print(molecule_area_c2h2)
        print(O_ss/rho_0)
        print((rho_0 - O_ss)/self.NAVA)
        print("dgffdgfdgfd")
        cross_section_area = (math.pi/4) * (d_outer**2 - d_inner**2) * rho_particle
        growth_rate_m_s = (carbon_atoms_flux * self.volume_carbon_atom) / cross_section_area
        print(growth_rate_m_s)
        return growth_rate_m_s

# Function to read Excel data
def read_excel_data(filename):
    """Read experimental data from Excel file"""
    try:
        # Read the Excel file with proper handling
        df = pd.read_excel(filename, header=1)  # Use row 2 as headers (0-indexed row 1)
        
        # Print column names to debug
        print("Column names:", df.columns.tolist())
        
        # Extract Liu et al. experimental data at 680°C
        liu_data = {
            'pressure': df.iloc[:, 0].dropna().values,  # Column A (first column)
            'growth_rate': df.iloc[:, 1].dropna().values * 1e-6  # Column B (second column), convert μm/s to m/s
        }
        
        # Extract Ma et al. simulation data at 680°C  
        ma_data = {
            'pressure': df.iloc[:, 3].dropna().values,  # Column D (fourth column)
            'growth_rate': df.iloc[:, 4].dropna().values * 1e-6  # Column E (fifth column), convert μm/s to m/s
        }
        
        print(f"Liu et al. data points: {len(liu_data['pressure'])}")
        print(f"Ma et al. data points: {len(ma_data['pressure'])}")
        print(f"Liu pressure range: {liu_data['pressure'].min():.0f} - {liu_data['pressure'].max():.0f} Pa")
        print(f"Ma pressure range: {ma_data['pressure'].min():.0f} - {ma_data['pressure'].max():.0f} Pa")
        print(f"Liu growth rate range: {liu_data['growth_rate'].min():.2e} - {liu_data['growth_rate'].max():.2e} m/s")
        print(f"Ma growth rate range: {ma_data['growth_rate'].min():.2e} - {ma_data['growth_rate'].max():.2e} m/s")
        
        return liu_data, ma_data
        
    except FileNotFoundError:
        print("Excel file not found.")
        return None, None
    
    except Exception as e:
        print(f"Error reading Excel file: {e}")
        import traceback
        traceback.print_exc()
        return None, None

# Main plotting function
def plot_pressure_vs_growth_rate():
    """Create plot of acetylene partial pressure vs growth rate"""
    
    # Temperature for the analysis (680°C = 953.15 K)
    temp_K = 680 + 273.15
    
   
    
   
    
   # Range of acetylene partial pressures (Pa) - from 2000 Pa to 10000 Pa
    pressure_range = np.linspace(2000, 10000, 20)  # 2,000 Pa to 10,000 Pa
    
    
    
    

    
    # Fixed C4H4 pressure
    P_FEEDSTOCK_c4h4_fixed = 0#456.645  # Pa
    
    
    
    
    
    
    # Particle diameter range for shaded area
    dp_range = np.linspace(5e-9, 20e-9, 15)  # 5-20 nm
    
    # Calculate growth rates for dp = 15 nm
    dp_main = 15e-9
    growth_rates_main = []
    
    print(f"Calculating growth rates for dp = {dp_main*1e9:.1f} nm at {temp_K-273.15:.0f}°C")
    ambient_temp=temp_K
    for P_c2h2 in pressure_range:
        
        gas_initial_composition={"C2H2" :   P_c2h2/101325, "C4H4" :   0,    "H2"   :   1e-35,"Ar"   :   0.900} ##molar fraction
        AA00=ch_sol.chemical_kinetics(gas=gas,flow_velocity=vel,chk_temp=ambient_temp,
                          reactor_pressure=pressure_atm,
                          inlet_composition=gas_initial_composition,
                          reactor_data={"length":reactor_length,
                                        "diameter":xdtube,
                                        "lengthstep":1000000,
                                        "reaction_mechanism":kinetics_mechanism_type
                                        })

        
        
        
        
        
        
        
        model = MaGrowthModel(dp=dp_main, P_FEEDSTOCK_c2h2=AA00[2]['C2H2']*101325, P_FEEDSTOCK_c4h4=AA00[2]['C4H4']*101325)
        growth_rate = model.calculate_growth_rate_per_step(temp_K)
        growth_rates_main.append(growth_rate)
    
    growth_rates_main = np.array(growth_rates_main)
    
    # Calculate growth rates for shaded area (different dp values)
    growth_rates_bounds = np.zeros((len(dp_range), len(pressure_range)))
    
    print("Calculating bounds for shaded area...")
    for i, dp in enumerate(dp_range):
        print(f"Processing dp = {dp*1e9:.1f} nm ({i+1}/{len(dp_range)})")
        for j, P_c2h2 in enumerate(pressure_range):
            gas_initial_composition={"C2H2" :   P_c2h2/101325, "C4H4" :   0,    "H2"   :   1e-35,"Ar"   :   0.900} ##molar fraction
            AA00=ch_sol.chemical_kinetics(gas=gas,flow_velocity=vel,chk_temp=ambient_temp,
                              reactor_pressure=pressure_atm,
                              inlet_composition=gas_initial_composition,
                              reactor_data={"length":reactor_length,
                                            "diameter":xdtube,
                                            "lengthstep": 1000000,
                                            "reaction_mechanism":kinetics_mechanism_type
                                            })
            model = MaGrowthModel(dp=dp, P_FEEDSTOCK_c2h2=AA00[2]['C2H2']*101325, P_FEEDSTOCK_c4h4=AA00[2]['C4H4']*101325)
            growth_rate = model.calculate_growth_rate_per_step(temp_K)
            growth_rates_bounds[i, j] = growth_rate
    
    # Get min and max bounds for shading
    growth_rates_min = np.min(growth_rates_bounds, axis=0)
    growth_rates_max = np.max(growth_rates_bounds, axis=0)
    
    # Read experimental data from Excel
    print("Reading experimental data from Excel...")
    liu_data, ma_data = read_excel_data('Ma_paper_digitized.xlsx')  # Updated filename
    
    # Create the plot
    plt.figure(figsize=(12, 12*6./8.18))
    plt.rcParams['font.family'] = 'Calibri'
    plt.rcParams['font.size'] = 24
    
    # Add shaded area for particle diameter range
    plt.fill_between(pressure_range, growth_rates_min, growth_rates_max,
                     color='#15C5C5', alpha=0.2, label=f'dp = 5-20 nm range')
    
    # Plot main curve for dp = 15 nm
    plt.semilogy(pressure_range, growth_rates_main, '--', linewidth=4,
               color='#2FAAE1', label=f'Ma model (dp = {dp_main*1e9:.0f} nm)', alpha=1)
    
    # Plot experimental data if available
    if liu_data is not None and ma_data is not None:
        plt.semilogy(liu_data['pressure'], liu_data['growth_rate'], 'o-', 
                   markersize=18, markerfacecolor='none', markeredgewidth=4,
                   color='#154357', label='Liu et al. exp (680°C)', alpha=1)
        
        plt.semilogy(ma_data['pressure'], ma_data['growth_rate'], 's--', 
                   markersize=18, markerfacecolor='none', markeredgewidth=4,
                   color='#724804', label='Ma et al. sim (680°C)', alpha=1)
    
    # Formatting
    plt.xlabel("Acetylene Partial Pressure [Pa]", fontsize=18, fontfamily='Calibri')
    plt.ylabel("Growth Rate [m/s]", fontsize=18, fontfamily='Calibri')
    
    # Set axis properties
    ax = plt.gca()
    ax.tick_params(axis='both', which='major', labelsize=18)
    ax.tick_params(axis='x', which='major', direction='out', bottom=True, top=False, width=2, length=8)
    ax.tick_params(axis='y', which='major', direction='out', left=True, right=False, width=2, length=8)
    ax.tick_params(axis='x', which='minor', direction='out', bottom=True, top=False, width=2, length=5)
    ax.tick_params(axis='y', which='minor', direction='out', left=True, right=False, width=2, length=5)
    
    # Spine styling
    ax.spines['bottom'].set_linewidth(2)
    ax.spines['left'].set_linewidth(2)
    ax.spines['top'].set_linewidth(2)
    ax.spines['right'].set_linewidth(2)
    
    # Grid and limits - adjusted for the specified pressure range
    plt.grid(False)
    plt.xlim(0, 11000)  # From 2000 Pa to 10000 Pa as requested
    plt.ylim(0, 9e-6)   # Adjusted to match experimental growth rates
    
    # Legend
    plt.legend(fontsize=16, loc='best')
    
    plt.tight_layout()
    
    # Define the save directory - using forward slashes and proper escaping
    save_dir = "C:/Users/rahba/Desktop/Phd/####Ph.D Thesis/####My Papers & Conferences&Proposals/####My Papers/4_Fourth_paper_CNT/##results/##Final_figures/##Paper_Final_Figures/surface_kinetics_verification/Ma_verifications/rate_vs_partial_pressure"
    
    # Alternative: Try to create directory step by step to debug
    try:
        # Create directory if it doesn't exist
        os.makedirs(save_dir, exist_ok=True)
        print(f"Directory created/verified: {save_dir}")
        
        # Test if directory is writable
        test_file = os.path.join(save_dir, "test_write.txt")
        with open(test_file, 'w') as f:
            f.write("test")
        os.remove(test_file)
        print("Directory is writable")
        
    except Exception as e:
        print(f"Error with directory: {e}")
        # Fallback to current directory
        save_dir = "."
        print("Saving to current directory instead")
    
    try:
        # Save the figure with full paths
        plt.savefig(os.path.join(save_dir, "cnt_growth_rate_vs_acetylene_pressure.jpg"), format='jpg', dpi=600, bbox_inches='tight')
        plt.savefig(os.path.join(save_dir, "cnt_growth_rate_vs_acetylene_pressure.svg"), format='svg', bbox_inches='tight')
        plt.savefig(os.path.join(save_dir, "cnt_growth_rate_vs_acetylene_pressure.tiff"), format='tiff', dpi=600, bbox_inches='tight')
        plt.savefig(os.path.join(save_dir, "cnt_growth_rate_vs_acetylene_pressure.pdf"), format='pdf', bbox_inches='tight')
        
        print(f"Figures successfully saved to: {save_dir}")
        
    except Exception as e:
        print(f"Error saving figures: {e}")
        # Try saving with simple filenames in current directory
        try:
            plt.savefig("cnt_growth_rate_vs_acetylene_pressure.jpg", format='jpg', dpi=600, bbox_inches='tight')
            plt.savefig("cnt_growth_rate_vs_acetylene_pressure.svg", format='svg', bbox_inches='tight')
            plt.savefig("cnt_growth_rate_vs_acetylene_pressure.tiff", format='tiff', dpi=600, bbox_inches='tight')
            plt.savefig("cnt_growth_rate_vs_acetylene_pressure.pdf", format='pdf', bbox_inches='tight')
            print("Figures saved to current working directory")
        except Exception as e2:
            print(f"Failed to save figures: {e2}")
    
    plt.show()
    
    return pressure_range, growth_rates_main, growth_rates_min, growth_rates_max

# Run the analysis
if __name__ == "__main__":
    pressure_range, growth_rates_main, growth_rates_min, growth_rates_max = plot_pressure_vs_growth_rate()
    
    print("\nAnalysis complete!")
    print(f"Pressure range: {pressure_range[0]:.0f} - {pressure_range[-1]:.0f} Pa")
    print(f"Growth rate range (dp=15nm): {growth_rates_main.min():.2e} - {growth_rates_main.max():.2e} m/s")

    print(f"Growth rate bounds: {growth_rates_min.min():.2e} - {growth_rates_max.max():.2e} m/s")
