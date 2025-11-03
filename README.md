# LongArithmetic

## Description

This project provides a `Large` class for performing arbitrary-precision integer arithmetic in C++.  
Large numbers are represented as arrays of digits in base \(10^9\).  

---

## Features

| Operation                     | Complexity |
|--------------------------------|-----------|
| Addition                       | O(n)      |
| Subtraction                    | O(n)      |
| Multiplication (small numbers) | O(n*m) or using `__int128` |
| Multiplication (large numbers) | O((n+m) log(n+m)) via FFT |
| Division                        | O(mult(n, m)) |
| Modulo                          | O(div(n, m)) |
| Exponentiation                  | O(log P * mult(n,n)) |
| Conversion to string            | O(n)      |
| Increment/Decrement             | O(n)      |
| Assignment                       | Move semantics |

---

## Constructors

The `Large` class provides three constructors:  
1. Default constructor — initializes to 0  
2. From `int64_t`  
3. From `std::string`  

---

## Usage

```cpp
#include "large-numbers.h"
#include <iostream>

int main() {
    Large a("123456789012345678901234567890");
    Large b(9876543210);

    Large result = a + b;
    std::cout << "Sum: " << result << "." << std::endl;

    return 0;
}
```