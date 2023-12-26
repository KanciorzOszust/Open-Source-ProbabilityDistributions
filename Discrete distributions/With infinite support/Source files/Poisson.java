public class Poisson extends MathLibrary{
    double lambda;
    Poisson(double lambda) {
        if (lambda <= 0) throw new IllegalArgumentException("lambda > 0");
        this.lambda = lambda;
    }

    double PMF(int k) {
        return (power(lambda, k) * euler(-1 + lambda)) / factorial(k);
    }

    double CDF(int k) {
        double value = euler(-1 * lambda);
        for (int i = 0; i < k; i++) {
            value += power(lambda, i) / factorial(i);
        }
        return value;
    }

    double Mean() {
        return lambda;
    }

    int Median() {
        return (int) (lambda + (1.0 / 3) - (1 / (50 * lambda)));
    }

    int Mode() {
        return (int) lambda;
    }

    double Variance() {
        return lambda;
    }

    double Skewness() {
        return 1 / root(lambda, 2);
    }

    double EXkurtosis() {
        return 1 / lambda;
    }

    double Entropy() {
        double wartosc = euler(-1 * lambda);
        for (int i = 0; i < 7; i++) {
            wartosc += (power(lambda, i) * log(factorial(i), 10)) / factorial(i);
        }
        return (lambda * (1 - log(lambda, 10))) + wartosc;
    }

    double MGF(double t) {
        return euler(lambda * (euler(t) - 1));
    }

    double PGF(double z) {
        return euler(lambda * (z - 1));
    }

    double FisherInformation() {
        return 1 / lambda;
    }
}
