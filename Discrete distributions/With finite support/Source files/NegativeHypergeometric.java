public class NegativeHypergeometric extends MathLibrary{
    int N;
    int K;
    int r;
    
    NegativeHypergeometric(int N, int K, int r) {
        if (N < 0 || K < 0 || r < 0) throw new IllegalArgumentException("N >= 0, K >= 0, n >= 0");
        if (K > N ) throw new IllegalArgumentException("K < N");
        if (r > (N - K)) throw new IllegalArgumentException("r < N - K");
        this.N = N;
        this.K = K;
        this.r = r;
    }

    double PMF(int k) {
        if (k < 0 || k > K) throw new IllegalArgumentException("0 <= k <= K");
        return (double) (binomial(k + r - 1, k) * binomial(N - r - k, K - k)) / binomial(N, K);
    }

    double Mean() {
        return r * ((double) K / (N - K + 1));
    }

    double Variance() {
        double part1 = ((double)(N + 1) * K) / ((N - K + 1) * (N - K + 2));
        double part2 = 1 - ((double) r / (N - K + 1));
        return r * part1 * part2;
    }
}
