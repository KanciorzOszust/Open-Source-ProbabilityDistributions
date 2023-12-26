public class ConwayMaxwellPoisson extends MathLibrary{
    double lambda;
    double nu;

    ConwayMaxwellPoisson(double lambda, double nu) {
        if (lambda <= 0) throw new IllegalArgumentException("lambda > 0");
        if (nu < 0) throw new IllegalArgumentException("nu >= 0");
        this.lambda = lambda;
        this.nu = nu;
    }

    private double NormalizationConstant(double a, double b) {
        double wartosc = 0;
        for (int i = 0; i < 7; i++) {
            wartosc += power(a, i) / doublePower(factorial(i), b);
        }
        return wartosc;
    }

    double PMF(int x) {
        return (power(lambda, x) / doublePower(factorial(x), nu)) * (1.0 / NormalizationConstant(lambda, nu));
    }

    double Mean() {
        double wartosc = 0;
        for (int i = 0; i < 7; i++) {
            double licznik = i * power(lambda, i);
            double mianownik = doublePower(factorial(i), nu) * NormalizationConstant(lambda, nu);
            wartosc += licznik / mianownik;
        }
        return wartosc;
    }

    double Variance() {
        double wartosc = 0;
        for (int i = 0; i < 7; i++) {
            double licznik = power(i, 2) * power(lambda, i);
            double mianownik = doublePower(factorial(i), nu) * NormalizationConstant(lambda, nu);
            wartosc += licznik / mianownik;
        }
        return wartosc - power(Mean(), 2);
    }

    double MGF(double t) {
        return NormalizationConstant(euler(t) * lambda, nu) / NormalizationConstant(lambda, nu);
    }

    double PGF(double t) {
        return NormalizationConstant(t * lambda, nu) / NormalizationConstant(lambda, nu);
    }
}
