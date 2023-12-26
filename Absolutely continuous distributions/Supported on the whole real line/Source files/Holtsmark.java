public class Holtsmark {
    double c;
    double mu;

    Holtsmark(double c, double mu) {
        if (c <= 0) throw new IllegalArgumentException("c > 0");
        this.c = c;
        this.mu = mu;
    }

    double Mean() {
        return mu;
    }

    double Median() {
        return mu;
    }

    double Mode() {
        return mu;
    }

    double Variance() {
        return Double.POSITIVE_INFINITY;
    }
}
