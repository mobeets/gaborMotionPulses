Y = spike rate
Yr = Y - Xw = neuron noise
X = stimulus (n-by-nw)

# Choice as a function of stimulus and noise

## 1. Choice = Y + noise

Choice is now driven to some extent both by stimulus and noise.

* CP(Y) contains both stimulus-driven and noise-driven changes in choice
* CP(Yr) contains only noise-driven changes in choice.

CP(Yr) <= CP(Y)

If CP(Yr) == CP(Y), then that's because the stimulus-driven effects cancel out and have no end contribution to choice.

## 2. Choice is dependent on Yr only.

* CP(Y) contains responses due to stimulus which are not choice related, so they only distract.

Thus CP(Y) <= CP(Yr).

# Choice with bimodal RF.

## 1. Choice = Y + noise

Suppose nw = 2.
wpos := w(w < 0) = 0
Ypos = Y - X*wpos
wpos := w(w > 0) = 0
Yneg = Y - X*wneg

Y = Ypos + Yneg + Yr

* As sum(w) increases, so does CP(Y)

This is because the noise in Y + noise is more dominated by Y.

* As sum(w) increases, CP(Yr) decreases

See above--stimulus drives choice more, so less effect from Yr.

* If sign(w[1]) == sign(w[2]) == 1, CP(Y - Yneg) == CP(Y)

Because Yneg == 0

* If sign(w[1]) == sign(w[2]) == 1, CP(Y - Ypos) == CP(Yr)

Because because Yr = Y - Ypos, since Yneg == 0.

* If prod(w) < 0, CP(Y) < CP(Y - Ypos) < CP(Yr)

Y - Ypos == Yr + Yneg, and both Yr and Yneg are correlated with choice.

* If prod(w) < 0, CP(Y) < CP(Y - Yneg) < CP(Yr)

Same as above basically.

* If prod(w) < 0 and abs(w[1]) > abs(w[2])

If sum(wpos) > sum(wneg), CP(Y - Ypos) < CP(Y - Yneg)

Because wpos drives choice more, subtracting off Ypos removes more choice information than removing Yneg.

* If CP(Y) < CP(Ypos + Yr), that must mean Yneg must be uncorrelated with choice.

Y = Ypos + Yneg + Yr

So if CP(Ypos + Yneg + Yr) < CP(Ypos + Yr), Yneg must be obscuring the discriminability of choice. 



Let C = Choice.
Suppose:
	A = C + eA*randn(ntrials,1)
	B = C + eB*randn(ntrials,1)
Then:
	CP(A+B)


## Requirements

Y = Ypos + Yneg + Yr

* CP(Yr) < CP(Y)
* CP(Ypos + Yr) > CP(Ypos + Yneg + Yr)
* corrcoef(Ypos, Yneg) << 0 (median around -0.75)

Ypos and Yneg are negatively correlated, so you want positive noise correlations.


Ehhhhhh

But given that CP(Ypos) > CP(Ypos + Yr), it seems like dominant correlation is from stimulus.
So adding in Yneg, i.e., CP(Ypos + Yneg + Yr), because Yneg is anti-correlated with YPos then this is going to lower overall CP to add it in, just because it's anti-correlated with the main driver of the CP.
