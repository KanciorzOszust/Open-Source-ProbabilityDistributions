public class HyperbolicSecant extends MathLibrary{
    double PDF(double x) {
        return 0.5 * (1 / hyperbolicCosine(PI / 2 * x));
    }

    double CDF(double x) {
        return (2 / PI) * arcTangent(euler(PI / 2 * x));
    }

    int Mean() {
        return 0;
    }

    int Median() {
        return 0;
    }

    int Mode() {
        return 0;
    }

    int Variance() {
        return 1;
    }

    int Skewness() {
        return 0;
    }

    int EXkurtosis() {
        return 2;
    }

    double Entropy() {
        return ln(4);
    }

    double MGF(double t) {
        if (absolute(t) < (PI / 2)) {
            return 1 / cosine(t);
        } else {
            return Double.NaN;
        }
    }
}
