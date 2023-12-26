public class Rademacher extends MathLibrary {
    double PMF(int k) {
        if (k == -1 || k == 1) return 0.5;
        return 0;
    }

    double CDF(double k) {
        if (k < - 1) return 0;
        else if (k >= 1) return 1;
        else if (k >= -1 && k < 1) return 0.5;
        else throw new IllegalArgumentException(" -1 <= k <= 1");
    }

    int Mean() {
        return 0;
    }

    int Median() {
        return 0;
    }

    int Variance() {
        return 1;
    }

    int Skewness() {
        return 0;
    }

    int EXkurtosis() {
        return -2;
    }

    double Entropy() {
        return ln(2);
    }

    double MGF(double t) {
        return hyperbolicCosine(t);
    }

    double CF(double t) {
        return cosine(t);
    }
}
