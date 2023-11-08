#!/usr/bin/python3
import numpy as np


N_STUDENTS = 29853
N_COLLEGES = 85
N_CHOICES = 5


def choice_one_time(n) -> list:
    choices = np.zeros(n, np.uint8)
    for i in range(N_CHOICES):
        rand_index = np.random.randint(1, N_COLLEGES)
        while rand_index in choices:
            rand_index = np.random.randint(1, N_COLLEGES)
        choices[i] = rand_index

    return choices


def make_data() -> list:
    array = []

    for _ in range(N_STUDENTS):
        array.append(choice_one_time(N_CHOICES))

    return array


def write_data(data: list):
    fh = open('./data/parents.txt', 'w')

    for s in data:
        for c in s:
            if c == s[-1]:
                print(c)
                fh.write(f"{c}\n")
                continue

            fh.write(f"{c},")
            print(c, end=',')


def main():
    write_data(make_data())


if __name__ == '__main__':
    main()
