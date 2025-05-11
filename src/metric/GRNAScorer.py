import csv
import ViennaRNA
import Bio.Seq

def calculate_gc_content(grna):
    gc_count = grna.count('G') + grna.count('C')
    return (gc_count / len(grna)) * 100

repeats = ['AAAAA', 'CCCCC', 'GGGG', 'UUUU']
def no_repetitive_bases(grna):
    return not any(rep in grna for rep in repeats)

def score_repetitive_bases(grna):
    if no_repetitive_bases(grna):
        return 4
    else:
        return -2.5*any(rep in grna for rep in repeats)

def uuu_not_in_seed(grna):
    seed_region = grna[-6:]
    return not 'UUU' in seed_region

def check_position_20(grna):
    return not grna[19] in ['C', 'U']

def evaluate_folding_energy(delta_g):
    return -1.9 >= delta_g >= -3.1

def score_folding_energy(delta_g):
    if evaluate_folding_energy(delta_g):
        return 4
    elif delta_g > -1.9:
        return (delta_g + 1.9) * -4
    else:
        return (delta_g + 3.1) * 4

def evaluate_duplex_stability(delta_g_duplex):
    return -20 >= delta_g_duplex >= -25

def score_duplex_stability(delta_g_duplex):
    if evaluate_duplex_stability:
        return 4
    elif delta_g_duplex > -20:
        return (delta_g_duplex + 20) * -0.5
    else:
        return (delta_g_duplex + 25) * 0.5

def assess_grna(grna):
    score = 0
    grna = grna.replace("T", "U")# RNA has U--but often stored in DB with T, mainly for checking repetitive bases
    #and UUU_Not_In_Seed region-----ViennaRNA treats T as U already
    gc_content = round(calculate_gc_content(grna), 2)
    has_optimal_gc = 30 <= gc_content <= 70
    score += 3 if has_optimal_gc else -3

    no_repeats = no_repetitive_bases(grna)
    score += score_repetitive_bases(grna)

    uuu_check = uuu_not_in_seed(grna)
    score += 1 if uuu_check else -1.5

    good_pos_20 = check_position_20(grna)
    score += 1 if good_pos_20 else -1.5

    # In this case folding energy is free energy
    structure, mfe = ViennaRNA.fold(grna)
    mfe_ok = evaluate_folding_energy(mfe)
    score += score_folding_energy(mfe)

    complement = Bio.Seq.reverse_complement(grna)
    duplex_energy = ViennaRNA.duplexfold(grna, complement).energy
    duplex_ok = evaluate_duplex_stability(duplex_energy)
    score += score_duplex_stability(duplex_energy)

    score = round(score, 2) #round to 2 decimal places
    likely_functional = (
            #uses duplex score but no sequences seem to have delta_g_duplex >= -22 and GC content >= 30
            score >= 9 and
            (has_optimal_gc and
            no_repeats and
            uuu_check and
            good_pos_20) and
            (mfe_ok or duplex_ok)
    )

    result = {
        'Score': score,
        'Likely_Functional': likely_functional,
        'gRNA': grna,
        'GC_Content': gc_content,
        'Has_Optimal_GC': has_optimal_gc,
        'No_Repetitive_Bases': no_repeats,
        'UUU_Not_In_Seed': uuu_check,
        'Good_Position_20': good_pos_20,
        'Folding_Energy': mfe,
        'Folding_Energy_OK': mfe_ok,
        'Duplex Energy': duplex_energy,
        'Duplex_Stability_OK': duplex_ok,
        'Structure': structure
    }

    return result

def batch_file(input_file, output_file):
    results = batch_assess_sort(input_file)
    with open(output_file, 'w', newline='') as csvfile:
        fieldnames = list(results[0].keys())
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for row in results:
            writer.writerow(row)

def batch_assess_sort(input_file):
    results = []
    with open(input_file, 'r') as file:
        for line in file:
            grna = line.strip()
            if grna:
                result = assess_grna(grna)
                results.append(result)
    results.sort(key=lambda x: x["Score"], reverse=True)
    return results

def batch_print(input_file):
    results = batch_assess_sort(input_file)
    for result in results:
        print_result_format(result)

def print_result_format(result):
    for key, value in result.items():
        if key == 'Likely_Functional':
            if value:
                print(f"{key}: \033[32mTrue\033[0m")  # Green True
            else:
                print(f"{key}: \033[31mFalse\033[0m") # Red False
        elif key == 'Score':
            if value >= 9:
                print(f"{key}: \033[32m{value}\033[0m")  # Green Score above 5
            else:
                print(f"{key}: \033[31m{value}\033[0m") # Red Score below 5
        else:
            print(f"{key}: {value}")

# Example batch usage:
# batch_assess('input_grnas.txt', 'output_metrics.csv')
def main():
    # Duplex_Stability_OK and GC_Content metrics make it hard to find a Likely_Functional gRNA
    batch_file('grnaSeqs.txt', 'grnaResults.csv')
    batch_print('grnaSeqs.txt')
    print("CHOPCHOP------"*9)
    batch_print('chopchop.txt')
    batch_file('chopchop.txt', 'chopchopResults.csv')


if __name__ == "__main__":
    main()
else:
    print("GRNA : Is intended to be executed and not imported.")
