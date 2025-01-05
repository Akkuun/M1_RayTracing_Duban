#ifndef __PROGRESSBAR_HPP
#define __PROGRESSBAR_HPP

#include <iostream>
#include <ostream>
#include <string>
#include <stdexcept>
#include <math.h>
#include <iomanip>

using namespace std;

class progressbar {

public:
    void update(double newProgress) {
        currentProgress = newProgress;
        amountOfFiller = (int) ((currentProgress / neededProgress) * (double) pBarLength);
    }
    string getRandomColor() {
        int colorCode;
        do {
            colorCode = rand() % 256; // Generate a random color code
        } while (colorCode == 0 || (colorCode >= 232 && colorCode <= 240)); // Exclude black and very dark colors
        return "\033[38;5;" + to_string(colorCode) + "m"; // ANSI 256-color foreground
    }



    void print(float remaningTime) {
        currUpdateVal %= pBarUpdater.length();

        // Calculate the amount of filler to display
        int amountOfFiller = static_cast<int>(pBarLength * currentProgress / neededProgress);

        // Apply random color
        cout << getRandomColor();

        // Bring cursor to start of line and print the progress bar
        cout << "\r" << firstPartOfpBar;
        for (int a = 0; a < amountOfFiller; a++) {
            cout << pBarFiller;
        }
        cout << pBarUpdater[currUpdateVal];
        for (int b = 0; b < (pBarLength - amountOfFiller); b++) {
            cout << " ";
        }

        cout << lastPartOfpBar << " "
              << static_cast<int>(100 * currentProgress / neededProgress) << "%";

        if (remaningTime != -1) { // If remaining time is provided
            // Format the remaining time to a fixed width
            cout << " Remaining Estimated time: " << std::setw(2) << std::fixed
                 << std::setprecision(1) << remaningTime << "s";
        }

        cout << "\033[0m" // Reset color
             << flush;

        currUpdateVal += 1;
    }


    std::string firstPartOfpBar = "[",
    lastPartOfpBar = "]",
            pBarFiller = "|",
            pBarUpdater = "/-\\|";

private:
    int amountOfFiller,
            pBarLength = 50,
    currUpdateVal = 0;
    double currentProgress = 0,
    neededProgress = 101;
};

#endif