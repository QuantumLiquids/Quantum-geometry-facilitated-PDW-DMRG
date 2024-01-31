import json
import matplotlib.pyplot as plt

# Define the initial values
Ly = 4
Lx = 32
ts = 1
td = -1
tsd_xy = 1
tsd_nn = 0
Uss = 3.8
Udd = 3.7
Usd = 4.0
Hole = 16

# Define the values of D to iterate over
D_values = [8000,10000]

for D in D_values:
    # Create the file path using f-string formatting
    file_path = f"../../data/nf{Ly}x{Lx}ts{ts}td{td}tsd_xy{tsd_xy}tsd_nn{tsd_nn}Uss{Uss}Udd{Udd}Usd{Usd}Hole{Hole}D{D}.json"

    # Load the data from the JSON file
    with open(file_path, 'r') as file:
        data = json.load(file)

    s_data = [entry for entry in data if entry[0][0] % (2*Ly) == 0]
    d_data = [entry for entry in data if entry[0][0] % (2*Ly) == 1]

    # Extract the x and y values
    x_values = [entry[0][0] // (2*Ly) for entry in s_data]  # Divide by Ly to get x
    n_s = [entry[1] for entry in s_data]
    n_d = [entry[1] for entry in d_data]

    plt.plot(x_values, n_s, marker='o', label=f'$n_s$, D = {D}')
    plt.plot(x_values, n_d, marker='x', label=f'$n_d$, D = {D}')

# Set the labels and title
plt.xlabel(r'$x$')
plt.ylabel(r'$\langle n(x)\rangle$')
plt.title('Charge Density Profile')

# Add a legend
plt.legend()

# Display the plot
plt.show()
