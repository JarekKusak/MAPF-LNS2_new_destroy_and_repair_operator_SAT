import matplotlib.pyplot as plt

def parse_sum_of_costs(logfile_path):
    sum_of_costs = []
    with open(logfile_path, "r", encoding="utf-8") as f:
        for line in f:
            if "[STAT] sum_of_costs after recomputation:" in line:
                parts = line.strip().split()
                # číslo bývá zpravidla na konci řádku
                cost_value = int(parts[-1])
                sum_of_costs.append(cost_value)
    return sum_of_costs

def plot_sum_of_costs(sum_of_costs):
    plt.plot(sum_of_costs)
    plt.xlabel("Iterace")
    plt.ylabel("Sum of cost")
    plt.title("Vývoj hodnoty Sum of cost v průběhu LNS iterací")
    plt.show()

if __name__ == "__main__":
    logfile = "log.txt" # soubor je ve stejném adresáři jako tento skript
    costs = parse_sum_of_costs(logfile)
    print("Načteno", len(costs), "iterací.")

    # vykreslení průběhu sum_of_costs
    plot_sum_of_costs(costs)