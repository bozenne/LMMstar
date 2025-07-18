* Version 1.2.0  (2025-07-25)
** User visible changes

** Internal change
- NEW FEATURE: ?? change in =effects= ??
- NEW FEATURE: a delta method via the =estimate= function can now be performed on the output when combining mixed model via =rbind= and =mlmm=. 
- NEW FEATURE: pooling has been added as a possible method, next to adjustment for multiple testing for the output of =anova=, =rbind=, and =mlmm=. 
- NEW FEATURE: a formula can now be input to =summarizeNA= in order to evaluate the missing data patterns per group.
- NEW FEATURE: an argument =force= can passed to =anova= to do a likelihood ratio test (LRT) without checking whether the models are nested nor whether they are fitted via REML.
- FIX: bug when having many variables in the formula of =summarize=.
- BREAK: =summarize= no longer evaluate correlations. This has been moved to a new function called =correlate=. 
#+BEGIN_SRC R :exports both :results output :session *R* :cache no
LMMstar.options(param.optimizer = c("init.cor" = "overall"))
#+END_SRC

** Internal change
- re-organize what the output of =anova= contain and corresponding methods (=coef=, =model.tables=, =vcov=, =iid=) to better align methods applied to the output of =anova=, =rbind=, and =mlmm=.
- BREAK: the optimization method of lmm is now checking whether the score w.r.t. to all parameters (instead of only the mean and variance parameters) is close to 0. This could impact the solution found by lmm.
- BREAK: initialization of the correlation parameters is now performed by averaging time-specific residual correlations instead of computing an overall correlation. This could impact the solution found by lmm. This behavior can be reverted using:
- NEW FEATURE: update method for the output of =rbind= and =mlmm=.


* Version 1.1.0  (2024-05-12)
** User visible changes
   - NEW FEATURE effects.lmm to evaluate average counterfactuals.
   - FIX: reparametrisation when using sigma() inside estimate
   - BREAK: update the calculation of the influence function with REML which may impact multiple testing adjustments. Use LMMstar.options(REML2ML=FALSE) to revert to old version.
   - CHANGE: Rename argument data into newdata in extractors such as logLik, information, model.frame, model.matrix, residuals, score, vcov
   - CHANGE: Rename argument keep.newdata into keep.data
   - CHANGE: Argument at for partial residual calculation (instead of using an attribute on argument type).
   - REMOVE: dependency on emmeans and sandwich.

** Internal change
   - partial derivative of the residuals w.r.t. model parameters can be exported as an attribute grad.
   - standard deviation for random effect computed via delta method instead of numerical differentiation


* Version 1.0.0 (2023-11-08)
** User visible changes
   - NEW FEATURE: support lme4 syntax for random intercepts.

** Internal change
   - substantial re-wrote of the functions for handling covariance patterns.
