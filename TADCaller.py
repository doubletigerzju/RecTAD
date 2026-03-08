import numpy as np
from scipy.signal import find_peaks
from scipy.stats import mannwhitneyu


class TADCaller:

    def __init__(self, min_window=10, max_window=30):
        self.min_window = min_window
        self.max_window = max_window


    def compute_score(self, M, window):

        N = M.shape[0]
        w = window
        h = w // 2
        score = np.zeros(N)

        for i in range(h, N - h):

            S = M[i-h:i+h, i-h:i+h]

            A = np.sum(S[:h, :h])
            B = np.sum(S[:h, h:])
            C = np.sum(S[h:, :h])
            D = np.sum(S[h:, h:])
            score[i] = (A + D) - (B + C)

        return score


    def detect_boundaries(self, score, distance):
        peaks, _ = find_peaks(score, distance=distance)
        return peaks


    def merge_peaks(self, peaks):

        if len(peaks) == 0:
            return peaks

        merged = []
        group = [peaks[0]]

        for p in peaks[1:]:

            if p == group[-1] + 1:
                group.append(p)
            else:
                merged.append(int(np.mean(group)))
                group = [p]

        merged.append(int(np.mean(group)))
        return np.array(merged)


    # -----------------------------
    # Statistical test for boundary
    # -----------------------------
    def test_boundary(self, M, i, window, alpha=0.05):

        N = M.shape[0]
        w = window
        h = w // 2

        if i - h < 0 or i + h >= N:
            return False

        S = M[i-h:i+h, i-h:i+h]

        A = S[:h, :h].flatten()
        B = S[:h, h:].flatten()
        C = S[h:, :h].flatten()
        D = S[h:, h:].flatten()

        AD = np.concatenate([A, D])
        BC = np.concatenate([B, C])

        effect = np.mean(AD) - np.mean(BC)

        if effect <= 0:
            return False

        stat, p = mannwhitneyu(AD, BC, alternative="greater")

        return p < alpha


    # -----------------------------
    # Filter boundaries by statistics
    # -----------------------------
    def filter_boundaries(self, M, boundaries):

        filtered = []
        w = self.max_window

        for b in boundaries:

            if self.test_boundary(M, b, w):
                filtered.append(b)

        return np.array(filtered)


    def multi_scale_peaks(self, M):

        N = M.shape[0]
        vote = np.zeros(N)
        window_count = 0

        for w in range(self.min_window, self.max_window + 1):

            score = self.compute_score(M, w)
            peaks = self.detect_boundaries(score, distance=w)

            for p in peaks:
                left = max(0, p - 1)
                right = min(N, p + 2)
                vote[left:right] += 1

            window_count += 1

        threshold = window_count / 2

        final_peaks = np.where(vote >= threshold)[0]
        final_peaks = self.merge_peaks(final_peaks)

        # statistical filtering
        final_peaks = self.filter_boundaries(M, final_peaks)

        return final_peaks


    def call_TADs(self, boundaries, N):

        TADs = []
        start = 0

        for b in boundaries:
            TADs.append((start, b))
            start = b

        TADs.append((start, N - 1))

        return TADs


    def fit(self, M):

        boundaries = self.multi_scale_peaks(M)
        TADs = self.call_TADs(boundaries, M.shape[0])

        return boundaries, TADs