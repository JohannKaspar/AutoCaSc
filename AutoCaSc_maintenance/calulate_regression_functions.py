import numpy
import matplotlib.pyplot as plt

def calculate_loeuf_coefficients():
    x = [0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5]
    y = [.99, 0.98, 0.95, 0.9, 0.7, 0.4, 0]

    mymodel = numpy.poly1d(numpy.polyfit(x, y, 4))

    myline = numpy.linspace(0, 0.5, 100)

    plt.scatter(x, y)
    plt.plot(myline, mymodel(myline))
    plt.xlabel("LOEUF")
    plt.ylabel("translated factor")
    plt.title("LOEUF translation to multiplication factor")
    plt.ylim(0, 1)
    plt.show()
    print(mymodel.coefficients)
    
def calculate_z_coefficients():
    x = [6, 5, 4, 3.09, 2.5, 2, 1, 0]
    y = [1.2, 1.2, 1.1, 1, 0.9, 0.8, 0.5, 0]
    mymodel = numpy.poly1d(numpy.polyfit(x, y, 3))
    myline = numpy.linspace(0, 6, 100)
    plt.scatter(x, y)
    plt.plot(myline, mymodel(myline))
    plt.xlabel("Z Score")
    plt.ylabel("translated factor")
    plt.title("Z Score translation to multiplication factor")
    plt.ylim(0, 1.1)
    plt.xlim(0, 3.1)
    plt.show()
    print(mymodel.coefficients)

def calculate_cadd_coefficients():
    x = [10, 20, 30, 40]
    y = [0.5, 0.9, 0.99, 1]

    mymodel = numpy.poly1d(numpy.polyfit(x, y, 3))

    myline = numpy.linspace(0, 60, 100)

    plt.scatter(x, y)
    plt.plot(myline, mymodel(myline))
    plt.xlabel("CADD Score")
    plt.ylabel("translated factor")
    plt.title("CADD translation to multiplication factor")
    plt.ylim(0, 1.1)
    plt.xlim(0, 60)
    plt.show()
    print(mymodel.coefficients)