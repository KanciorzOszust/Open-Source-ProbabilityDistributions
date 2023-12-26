public class Degenerate extends MathLibrary {
    double k0;

    Degenerate(double k0) {
        this.k0 = k0;
    }

    int PMF(double x) {
        if (k0 == x) return 1;
        else return 0;
    }

    int CDF(double x) {
        if (x < k0) return 0;
        else return 1;
    }

    double Mean() {
        return k0;
    }

    double Median() {
        return k0;
    }

    double Mode() {
        return k0;
    }

    int Variance() {
        return 0;
    }

    int Entropy() {
        return 0;
    }

    double MGF(double t) {
        return euler(k0 * t);
    }
}
