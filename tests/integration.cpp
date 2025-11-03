#include <iostream>
#include <large-numbers.h>
#include <conio.h>

int main()
{
    setlocale(LC_ALL, "RU");

    Large first_number;
    Large second_number;

    std::cout << "Enter your first number: "; std::cin >> first_number;
    std::cout << "Enter your second number: "; std::cin >> second_number;

    std::cout << std::endl;
    std::cout << "The result of `first + second` = " << first_number + second_number;

    std::cout << std::endl;
    std::cout << "The result of `first - second` = " << first_number - second_number;

    std::cout << std::endl;
    std::cout << "The result of `first * second` = " << first_number * second_number;

    std::cout << std::endl;
    std::cout << "The result of `first / second` = " << first_number / second_number;

    std::cout << std::endl;
    std::cout << "The result of `first % second` = " << first_number % second_number;

    std::cout << std::endl;
    std::cout << "The result of `first ^ second` = " << pow(first_number, second_number);

    _getch();
}