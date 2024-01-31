import json
import matplotlib.pyplot as plt

# Define the initial values
Ly = 3
Lx = 24
ts = 1
td = -1
tsd_xy = 1
tsd_nn = 0
Uss = 3.8
Udd = 3.7
Usd = 4.0
Hole = 8
D_values = [10000,12000,14000,16000]

for D in D_values:
    # Create the file path using f-string formatting
    file_path = f"../data/zzsf{Ly}x{Lx}ts{ts}td{td}tsd_xy{tsd_xy}tsd_nn{tsd_nn}Uss{Uss}Udd{Udd}Usd{Usd}Hole{Hole}D{D}.json"
    
    # Load the data from the JSON file
    with open(file_path, 'r') as file:
        data = json.load(file)
    
    # Filter the data based on entry[0][0] == Lx * Ly / 2
    filtered_data = [entry for entry in data if (entry[0][0] == Lx * Ly / 2 and (entry[0][1] - entry[0][0]) % (2*Ly) == 0)]
    
    # Extract the x and y values
    x_values = [(entry[0][1] - entry[0][0]) // (2*Ly) for entry in filtered_data]  # Divide by Ly to get x
    y_values = [ -entry[1] for entry in filtered_data]
    
    # Plot the data on a logarithmic scale
    #y_values = [((-1) ** x) * y for x, y in zip(x_values, y_values)]
    plt.semilogy(x_values, y_values, marker='o')


plt.xlabel(r'$x$')
plt.ylabel(r'-$\langle S^z(0) S^z(x) \rangle$')
#plt.ylim(1e-5, 1e-1)  
plt.show()
