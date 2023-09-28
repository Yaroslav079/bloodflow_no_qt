import os
import csv
import matplotlib.pyplot as plt

menu_options = {
    1: 'Write errors from filter log to csv file. After number of option type path to log file',
    2: 'Create plot of filter errors. After number of option type path to error csv file',
    3: 'Write errors from bayes log to csv file. After number of option type path to log file',
    4: 'Create plot of bayes errors. After number of option type path to error csv file',
    42: 'Exit',
}

def print_menu():
    for key in menu_options.keys():
        print (key, '--', menu_options[key])

def errors_from_log(log_path, pattern, part_of_error_in_pattern):
    errors = []
    with open(log_path, 'r') as log_file:
        lines = log_file.readlines()
        for line in lines:
            if line.find(pattern) != -1:
                errors.append(line.split(pattern)[part_of_error_in_pattern])
    return errors

def plot_error(path, steps_in_iteration):
    errors = []
    with open(path, 'r') as csvfile:
        data = csv.reader(csvfile, delimiter = ',')
        for row in data:
            errors.append(float(row[0]) * 100.0)
    iter_nums = [*range(1, steps_in_iteration * len(errors) + 1, steps_in_iteration)]
    plt.plot(iter_nums, errors, 'bo-')
    plt.grid()
    plt.xlabel('Number of iteration')
    plt.ylabel('Error in percent')
    # plt.yscale('log')
    plt.show()

def list_to_csv(list, path):
    with open(path, 'w') as file:
        file.writelines(list)

if __name__ == '__main__':
    while True:

        print_menu()
        user_input = input('Choose option : ')
        input_words = user_input.split(' ')
        option = (int)(input_words[0])

        if option == 1:
            if os.path.exists(input_words[1]) == False:
                print('Choose existing log file')
            else:
                path_to_csv = input_words[1].split('.')[0] + '.csv'
                list_to_csv(errors_from_log(input_words[1], 'norm : ', 1), path_to_csv)

        elif option == 2:
            if os.path.exists(input_words[1]) == False:
                print('Choose existing log file')
            else:
                plot_error(input_words[1], 8)

        elif option == 3:
            if os.path.exists(input_words[1]) == False:
                print('Choose existing log file')
            else:
                path_to_csv = input_words[1].split('.')[0] + '.csv'
                list_to_csv(errors_from_log(input_words[1], 'Target registered : ', 1), path_to_csv)

        elif option == 4:
            if os.path.exists(input_words[1]) == False:
                print('Choose existing log file')
            else:
                plot_error(input_words[1], 1)

        elif option == 42:
            exit()

        else:
            print('Invalid option. Choose proper one')
