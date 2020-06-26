# Generalized Additive Model Selection

a fork by [EDG](https://egenn.github.io)  
***work in progress***

Reference: [Chouldechova A and Hastie T, 2015](https://arxiv.org/abs/1506.03850)

* Data-derived defaults for `degrees` and `dfs`
* New arguments `min.degree`, `max.degree`, `min.df`, `max.df`
* `failsafe` option that fits a `glm` instead (will be changed to `glmnet`).

The above try to prevent the algorithm from hanging or failing, which happens with bad values for `degrees` or `dfs`. These little fixes should probably be improved. They were implemented to allow running gamsel within another algorithm on a large number of experiments.
