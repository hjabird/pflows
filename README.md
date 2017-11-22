# pflows
Numerical implementations of papers.

This is not a high quality unit tested library or set of reference implementations or any such like.
It is a way of me keeping my potential flow theory based code organised.

* Sclavounos 1987 - "An unsteady lifting line theory" - Implementation (currently not working)
* Guermond & Sellier 1991 - "A unified unsteady lifting line theory" - Some equations

## Prerequisits:

Libraries:
* HBTK - See github.com/hjabird/HBTK - integration & GNUPlot integration
* Eigen3 - Google it - linear algebra
* Boost - special functions (Using an interface, so this could be changed without much pain).
* GNUPlot - plotting program

THis has only been built using MSVC from Visual Studio. You'll need a C++17 compatible compiler for everything to work (probably). 
If you want to work on this, adding the Eigen's debugging NATVIS to your Visual Studio installation is useful. Ask the Google.
