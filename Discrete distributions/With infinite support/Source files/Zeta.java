public class Zeta extends MathLibrary{
    double s;

    Zeta(double s) {
        if (s <= 1) throw new IllegalArgumentException("s > 1");
        this.s = s;
    }

    double PMF(int k) {
        return 1.0 / doublePower(k, s) / zeta(s);
    }

    double CDF(int k) {
        return (double) generalizedHarmonicNumber(k, s) / zeta(s);
    }

    double Mean() {
        if (s > 2) {
            return (double) zeta(s - 1) / zeta(s);
        } else {
            return 0;
        }
    }

    int Mode() {
        return 1;
    }

    double Variance() {
        if (s > 3) {
            return (double) (zeta(s) * zeta(s - 2) - power(zeta(s - 1), 2)) / power(zeta(s), 2);
        } else {
            return 0;
        }
    }

    double Entropy() {
        double wartosc = 0;
        for (int i = 1; i < 7; i++) {
            wartosc += (1.0 / doublePower(i, s) / zeta(s)) * log(doublePower(i, s) * zeta(s), 10);
        }
        return wartosc;
    }
}
