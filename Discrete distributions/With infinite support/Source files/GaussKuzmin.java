public class GaussKuzmin extends MathLibrary {
    double PMF(int k) {
        double x = 1 - (1.0 / power(k + 1, 2));
        return -1 * log(x, 2);
    }

    double CDF(int k) {
        return 1 - log((k + 2) / (k + 1), 2);
    }

    double Mean() {
        return Double.POSITIVE_INFINITY;
    }

    int Median() {
        return 2;
    }

    int Mode() {
        return 1;
    }

    double Variance() {
        return Double.POSITIVE_INFINITY;
    }

    double Entropy() {
        return 3.432527514776;
    }
}
