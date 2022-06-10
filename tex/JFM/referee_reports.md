# Referee 1 #

This paper examines the classic problem of thermal convection in a cylindrical annulus with sloping end walls (as if extracted from a sphere).  The setup uses Busse’s formulation where the 3D annulus is  approximated by suitable Cartesian 2D equations involving a beta term to represent the vortex stretching induced by the sloping end walls.  The main point of this paper to examine the effectiveness of the CE2 cumulant expansion (or truncation) in Direct Statistical Simulations (DSS) for this problem.  This problem has multiple basins of attraction even at moderate nonlinearity and therefore has interestingly complex dynamics, such as relaxation oscillations between jet-dominated behavior and small-scale turbulence dominance.   This paper concludes that the CE2 DSS *can* obtain results in agreement with (much more computationally expensive) direct numerical simulations (DNS) but it is necessary to be very careful about the information built into the initial conditions provided.  The authors bring up the concepts of “maximal knowledge” initial conditions (from DNS), that have the correct symmetries built in, and “maximal ignorance” initial conditions that are random and have no preordained information (and in-between initial conditions that may be random but with some symmetry bias).  Whether the CE2 DSS successfully mimics the DNS results is sensitive to these initial conditions.

This paper is relevant, clearly written and worth adding to the literature.  It is an interesting exposé of the foibles and difficulties associated with attempting to do less computation in highly nonlinear problems.  I certainly think it deserves to be published, pretty much as is. 

My only major comment comments are: 

1.  This paper does not make its place in the hierarchy of the papers by (combinations of) the authors at all clear.  Readers may be mystified as to why CE2 is being discussed when Marston et al (2018) and Tobias et al (2018) have already shown that Generalized QuasiLinear and CE2.5 or CE3 (other levels of DSS) are required to give accurate results.   The paper here *does* provide new insight on the role of initial conditions and symmetries for CE2 that is not in the other papers, but this is not made clear in the title, abstract or introduction.   I think rectification of this omission would really help the reader who has some acquaintance with the problem.

2.  Do the authors in conclusion wish to offer guidance as to how and when to use CE2 then?  If one is starting a new problem where not much is known, how should a researcher proceed?  Some assessment of when and where CE2 is likely to be useful would be helpful and again, justify the place of this paper in amongst the others a bit more.   (However, I do appreciate the discussion of how to *improve* beyond CE2 in the Discussion).

The rest of my comments are very minor, but here they are, more in chronological order than any order of importance:

Intro:
====

I would use “Buses annulus” rather than “Buses Annulus” as the authors did in previous papers.

**response: fixed.**

The paper uses an abundance of semicolons in the “…clause; clause …”.  There are numerous per paragraph and sometimes even nested versions.  Whilst I appreciate that this is not incorrect grammatically, I find that it does not help the readability having so many.

“CE2 has been shown to be effective …tightly coupled … shear flow”.  Needs a reference.
**response: references added**

Section 2:
=======

The Cartesian directions are not perfectly clear.   The domain is 2D in (x,y) and Omega is in the third (z) direction.

**This paragraph was re-written to make clear the domain**

“The system is governed by four dimensionless parameters …”.   These should be defined here and not pushed to a reference.

**Descriptions of each parameter have been added**

“We have found fully explicit time stepping to be more stable for CE2”.   This odd but perhaps explanation is beyond the desired scope of this paper.

“… we found it necessary to remove negative eigenvalues …”.   This is also a bit odd for this uninitiated in DSS.   The two last statements here in this report give the reader a bit of a sense of unease.  Maybe some explanation is required for anxiety reduction?

**We have added an additional section describing the rank instability of CE2 and its generation of negative eigenvalues**

Section 3:
========

“which also produces a 3-jet solution; an unsurprising result given that CE2 is fundamentally a quasilinear
theory”.  This semi-colon really should be a comma!

**Fixed**

“if the symmetry of the initial condition is sffciently biased.”    Some indication of the value of lambda that is sufficient would be helpful.

“This type of behaviour has been described in terms of predator-prey dynamics with the convective turbulence taking the role of the prey and the zonal  flows acting as a predator.”  Needs reference.

**citation added**

“Remarkably, CE2 is able to continue relaxation oscillations from this state; we believe that this is the rst time a quasilinear model has reproduced such complicated behaviour; previously they have been …”.   Two semi-colons in one sentence!

**this sentence has been replaced to address concerns of other referees**

“This solution clearly does not replicate the correct number of jets — though it does show some indication of bursting behaviour, as the strength of each zonal jet waxes and wanes in response to the driving.”  Might this solution eventually relax to the correct number of jets?   I’m surprised that it has too many, as more jets are harder to obtain than less jets.

**We have no evidence that the CE2 solution will relax further. Given that the CE2 equations are stiffer than DNS, the timescales over which significant changes in evolution are much shorter for the former. We are not quite certain what the referee is referring to with regard to the fact that more jets are harder to obtain than fewer.**


# Referee 2 #

Comments to the Author
The Busse annulus model of rotating convection has a rich dynamics, with zonal jets being formed, and these can give rise to relaxation oscillations as the jets suppress convection, despite the jets being created by the convection. A further complication is that the number of zonal jets formed depends on the initial conditions, because a number of different stable zonal flow configurations can exist. These phenomena are usually investigated by direct numerical simulation (DNS), but in this paper a relatively new technique, direct statistical simulation (DSS) is also used. The main objective is to see whether DSS can capture the essential behaviour found in the more time-consuming DNS approach. A secondary objective is to use the comparison between the DSS results and the DNS to identify what are the crucial features necessary to explain the zonal flow behaviour, in particular whether there are types of interactions which are not crucial. If this is the case, it opens up a simplified, and potentially less computationally expensive way of investigating these fundamental types of fluid motion. The results show that DSS does a pretty good job of reproducing the main features of Busse annulus behaviour, indicating that the omitted eddy-eddy interactions in the eddy equations don’t play a key role in this problem, but nevertheless do affect which basin of attraction you end up in.
<p>
While the manuscript is generally clearly written, I have a few points, mostly minor, which I think will improve the presentation.
<p>
While I understand the authors want to keep the paper as short as possible, I think it would be useful to write down the equations for c_{psi zeta} and c_{psi theta} as these are required to solve (2.7) and (2.8). I expect the required expressions can be derived from c_{zeta zeta} and c_{zeta theta} , but not all readers will be familiar with DSS, so it would help readers to understand what has been done if these formulae were stated explicitly. Also, a little more detail about the symmetries would be useful: presumably c_{theta zeta} = c_{zeta theta}, but why not say so explicitly?
<p>

**We added some additional text explaining these issues.**
Page 2 line 6 up:  better would be ‘y-direction (d), so 0 <= y,=1, a typical timescale …  . You need to specify that the boundaries are at y=0 and y=1 somewhere.
<p>

**fixed**
Page 5 line 6: ‘an unsurprising result … quasilinear theory.’ I think this remark comes from the fact that the Reynolds stress derived from the linear theory perturbations is symmetric about y=0.5, which the 3 jet solution is, but the two jet solution is not. However, this is not mentioned in the text and I am not sure all readers will be aware of this, so I think this should be mentioned here to increase clarity.
<p>
**here, we simply meant that since both DSS and the single-realization QL from Tobias et al (2018) are both quasi-linear theories, it is not surprising that they both fail in the same way. We have clarified this point in the text**

Page 6 equation 3.1. Not sure why you have a subscript 0 for y here, and it would be better to write pi y / L_y for both terms.
<p>

**Fixed**
Page 6 line 2 up. ‘settling down into a five jet solution’. Do you know whether DNS always settles down to the 5 jet solution irrespective of initial conditions? The paragraph is written as though there is only one final solution in DNS for a given set of parameters, but is it possible there can be several different jet solutions in DNS if different initial conditions are tried? Maybe you don’t know the answer to this, but a comment on it might be helpful.
<p>

**We have only ever found five jet solutions; this is consistent with similar parameters in Rotvig and Jones (JFM 2006). It is certainly not out of the question that another solution may be stable, however.**
Page 8, Figure 5. The orange and blue curves overlap in the top and middle panel, so orange is invisible. Note the overlap in the figure caption.
**Noted.**


# Referee 3 #

General comment 1
**We have added references to both S3T and the work of Thompson (1970). For the latter, we point out that the Busse annulus is a bit different from Thompson's shearing convection state. In the Busse annulus system, the symmetry is broken by the $\beta$ parameter. It does not feature sponaneous symmetry breaking pitchfork bifurcation to the shearing state.**

General comment 2
**Steve, do you want to have a look at this?**

General comment 3
**While this is an intriguing idea, and one that should absolutely be pursued, we do not believe we can address it in this work due to the space constraints.**

Specific comments

1. **we have added a reference to Farrell and Ioannou (2019) for the effect of jet-turbulence interaction at $\beta = 0$.

2. **References added, as noted above**

3. **References added**

4. 

5. sec2 **We have rewritten the description to aid understanding.**

6. eq. 2.6 **That is correct; we have corrected the equation.**

7. "this implies..." **This refers to the boundary conditions, not the midplane. We have clarified the point.**

8. **We have addressed the issue of removing negative eigenvalues in response to another referee's comments.**


9. "completely uncorrelated" **We have clarified this initial condition.** 

10. "Addtional question" **That is correct. We do not add any additional noise; the CE2 system chooses between the two via the rotation.**

11. "typo" **Fixed.**

12. "Comments on table 1" **We have corrected the typos. The physical meaning of runs with C=0 is that the Ekman number becomes asymptotically small such that the corresponding Ekman layer friction becomes unimportant.**

13. "Started from maximum ignorance..." **Not shown added**

14. "Formatting" **Fixed**

15. "I think the authors..." **That is correct. We have clarified in the text.**

16. "The authors don't comment..." **We have simplified this and and clarified.**

17. "It is interesting to me..." **STEVE/BRAD/KEATON, this is actually interesting. Need to think on it a bit. Given our speculation that EENL just adds dissipation, why is this not more energetic? It may be a non-obvious interaction with the convection as the referee points out, but I'm not sure.

18. "The discussion surrounding..." **STEVE/BRAD/KEATON, it is surprising. Do we want to try to provide more info for some kind of physical interpretation, as this referee suggests?

19. "Is there a good reason..." **We provided the enstrophy spectrum because it connects with the dynamical variable of interest $c_{\zeta \zeta}$. 
