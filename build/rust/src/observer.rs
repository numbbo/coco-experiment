use coco_sys::coco_observer_t;
use std::ffi::{CStr, CString};

/// An observer to log results in COCO's data format.
///
/// Can be provided to [Suite::next_problem] and it will
/// automatically be attached to the returned problem.
pub struct Observer {
    pub(crate) inner: *mut coco_observer_t,
}

pub enum ObserverName {
    Bbob,
    BbobBiobj,
    Toy,
    None,
}
impl ObserverName {
    fn as_str(&self) -> &'static str {
        match self {
            ObserverName::Bbob => "bbob",
            ObserverName::BbobBiobj => "bbob-biobj",
            ObserverName::Toy => "toy",
            ObserverName::None => "no-observer",
        }
    }
}

impl Observer {
    /// Creates a new observer.
    ///
    /// # observer_options
    /// A string of pairs "key: value" used to pass the options to the observer. Some
    /// observer options are general, while others are specific to some observers. Here we list only the general
    /// options, see observer_bbob, observer_biobj and observer_toy for options of the specific observers.
    /// - "result_folder: NAME" determines the folder within the "exdata" folder into which the results will be
    /// output. If the folder with the given name already exists, first NAME_001 will be tried, then NAME_002 and
    /// so on. The default value is "default".
    /// - "algorithm_name: NAME", where NAME is a short name of the algorithm that will be used in plots (no
    /// spaces are allowed). The default value is "ALG".
    /// - "algorithm_info: STRING" stores the description of the algorithm. If it contains spaces, it must be
    /// surrounded by double quotes. The default value is "" (no description).
    /// - "number_target_triggers: VALUE" defines the number of targets between each 10^i and 10^(i+1)
    /// (equally spaced in the logarithmic scale) that trigger logging. The default value is 100.
    /// - "target_precision: VALUE" defines the precision used for targets (there are no targets for
    /// abs(values) < target_precision). The default value is 1e-8.
    /// - "number_evaluation_triggers: VALUE" defines the number of evaluations to be logged between each 10^i
    /// and 10^(i+1). The default value is 20.
    /// - "base_evaluation_triggers: VALUES" defines the base evaluations used to produce an additional
    /// evaluation-based logging. The numbers of evaluations that trigger logging are every
    /// base_evaluation * dimension * (10^i). For example, if base_evaluation_triggers = "1,2,5", the logger will
    /// be triggered by evaluations dim*1, dim*2, dim*5, 10*dim*1, 10*dim*2, 10*dim*5, 100*dim*1, 100*dim*2,
    /// 100*dim*5, ... The default value is "1,2,5".
    /// - "precision_x: VALUE" defines the precision used when outputting variables and corresponds to the number
    /// of digits to be printed after the decimal point. The default value is 8.
    /// - "precision_f: VALUE" defines the precision used when outputting f values and corresponds to the number of
    /// digits to be printed after the decimal point. The default value is 15.
    /// - "precision_g: VALUE" defines the precision used when outputting constraints and corresponds to the number
    /// of digits to be printed after the decimal point. The default value is 3.
    /// - "log_discrete_as_int: VALUE" determines whether the values of integer variables (in mixed-integer problems)
    /// are logged as integers (1) or not (0 - in this case they are logged as doubles). The default value is 0.
    pub fn new(name: ObserverName, options: &str) -> Option<Observer> {
        let name = CString::new(name.as_str()).unwrap();
        let options = CString::new(options).unwrap();

        let inner = unsafe { coco_sys::coco_observer(name.as_ptr(), options.as_ptr()) };

        if inner.is_null() {
            None
        } else {
            Some(Observer { inner })
        }
    }

    /// Prints where the result is written to.
    pub fn result_folder(&self) -> &str {
        unsafe {
            CStr::from_ptr(coco_sys::coco_observer_get_result_folder(self.inner))
                .to_str()
                .unwrap()
        }
    }
}

impl Drop for Observer {
    fn drop(&mut self) {
        unsafe {
            coco_sys::coco_observer_free(self.inner);
        }
    }
}
