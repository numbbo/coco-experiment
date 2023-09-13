use coco_sys::coco_problem_t;
use std::{ffi::CStr, marker::PhantomData, ops::RangeInclusive};

use crate::suite::Suite;

/// A specific problem instance.
///
/// Instances can be optained using [Suite::next_problem]
/// and [Suite::problem_by_function_dimension_instance].
pub struct Problem<'suite> {
    pub(crate) inner: *mut coco_problem_t,
    pub(crate) function_idx: usize,
    pub(crate) instance_idx: usize,
    pub(crate) dimension_idx: usize,
    pub(crate) _phantom: PhantomData<&'suite Suite>,
}

unsafe impl Send for Problem<'_> {}

impl Problem<'_> {
    /// Returns the ID of the problem.
    ///
    /// For the `toy` suite this is
    /// - `{function-name}_d{dimension}`
    ///
    /// For `bbob` it is
    /// - bbob_f{function-index}_i{instance}_d{dimension}
    pub fn id(&self) -> &str {
        unsafe {
            CStr::from_ptr(coco_sys::coco_problem_get_id(self.inner))
                .to_str()
                .unwrap()
        }
    }

    /// Returns the name of the problem.
    pub fn name(&self) -> &str {
        unsafe {
            CStr::from_ptr(coco_sys::coco_problem_get_name(self.inner))
                .to_str()
                .unwrap()
        }
    }

    /// Returns the index of the problem.
    pub fn function_index(&self) -> usize {
        self.function_idx
    }

    /// Returns the dimension index of the problem.
    pub fn dimension_index(&self) -> usize {
        self.dimension_idx
    }

    /// Returns the instance of the problem.
    pub fn instance_index(&self) -> usize {
        self.instance_idx
    }

    /// Evaluates the problem at `x` and returns the result in `y`.
    ///
    /// The length of `x` must match [Problem::dimension] and the
    /// length of `y` must match [Problem::number_of_objectives].
    pub fn evaluate_function(&mut self, x: &[f64], y: &mut [f64]) {
        assert_eq!(self.dimension(), x.len());
        assert_eq!(self.number_of_objectives(), y.len());

        unsafe {
            coco_sys::coco_evaluate_function(self.inner, x.as_ptr(), y.as_mut_ptr());
        }
    }

    /// Evaluates the problem constraints in point x and save the result in y.
    ///
    /// The length of `x` must match [Problem::dimension] and the
    /// length of `y` must match [Problem::number_of_constraints].
    pub fn evaluate_constraint(&mut self, x: &[f64], y: &mut [f64]) {
        assert_eq!(self.dimension(), x.len());
        assert_eq!(self.number_of_constraints(), y.len());

        unsafe {
            coco_sys::coco_evaluate_constraint(self.inner, x.as_ptr(), y.as_mut_ptr());
        }
    }

    /// Returns true if a previous evaluation hit the target value.
    pub fn final_target_hit(&self) -> bool {
        unsafe { coco_sys::coco_problem_final_target_hit(self.inner) == 1 }
    }

    /// Returns the optimal function value + delta of the problem
    pub fn final_target_value(&self) -> f64 {
        unsafe { coco_sys::coco_problem_get_final_target_fvalue1(self.inner) }
    }

    /// Returns the optimal function value of the problem
    ///
    /// To check whether the target has been reached use [[Problem::final_target_value]]
    /// or [[Problem::final_target_hit]] instead.
    pub fn best_value(&self) -> f64 {
        unsafe { coco_sys::coco_problem_get_best_value(self.inner) }
    }

    /// Returns the dimension of the problem.
    pub fn dimension(&self) -> usize {
        unsafe {
            coco_sys::coco_problem_get_dimension(self.inner)
                .try_into()
                .unwrap()
        }
    }

    /// Returns the number of objectives of the problem.
    pub fn number_of_objectives(&self) -> usize {
        unsafe {
            coco_sys::coco_problem_get_number_of_objectives(self.inner)
                .try_into()
                .unwrap()
        }
    }

    /// Returns the number of constraints of the problem.
    pub fn number_of_constraints(&self) -> usize {
        unsafe {
            coco_sys::coco_problem_get_number_of_constraints(self.inner)
                .try_into()
                .unwrap()
        }
    }

    /// Returns the numver of integer variables of the problem.
    ///
    /// The first `n` variables will be integers then.
    /// Returns `0` if all variables are continuous.
    pub fn number_of_integer_variables(&self) -> usize {
        unsafe {
            coco_sys::coco_problem_get_number_of_integer_variables(self.inner)
                .try_into()
                .unwrap()
        }
    }

    /// Returns the upper and lover bounds of the problem.
    pub fn get_ranges_of_interest(&self) -> Vec<RangeInclusive<f64>> {
        let dimension = self.dimension() as isize;
        unsafe {
            let smallest = coco_sys::coco_problem_get_smallest_values_of_interest(self.inner);
            let largest = coco_sys::coco_problem_get_largest_values_of_interest(self.inner);

            (0..dimension)
                .into_iter()
                .map(|i| (*smallest.offset(i))..=(*largest.offset(i)))
                .collect()
        }
    }

    /// Returns how often this instance has been evaluated.
    pub fn evaluations(&self) -> u64 {
        unsafe {
            #[allow(clippy::useless_conversion)]
            coco_sys::coco_problem_get_evaluations(self.inner)
                .try_into()
                .unwrap()
        }
    }

    /// Returns how often this instances constrants have been evaluated.
    pub fn evaluations_constraints(&self) -> u64 {
        unsafe {
            #[allow(clippy::useless_conversion)]
            coco_sys::coco_problem_get_evaluations_constraints(self.inner)
                .try_into()
                .unwrap()
        }
    }

    /// Writes a feasible initial solution into `x`.
    ///
    /// If the problem does not provide a specific solution,
    /// it will be the center of the problem's region of interest.
    pub fn initial_solution(&self, x: &mut [f64]) {
        assert_eq!(self.dimension(), x.len());
        unsafe {
            coco_sys::coco_problem_get_initial_solution(self.inner, x.as_mut_ptr());
        }
    }
}

impl Drop for Problem<'_> {
    fn drop(&mut self) {
        unsafe {
            coco_sys::coco_problem_free(self.inner);
        }
    }
}
