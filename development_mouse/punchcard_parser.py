from typing import *
import yaml
import os
import luigi
import cytograph as cg
import development_mouse as dm
import logging
from collections import defaultdict
import copy

# those are the analyses allowed, if a kind of analysis is not here cannot be run using the analysis submodule
require_type_dict = {"Level1": dm.Level1, "PerformAnalysis": dm.PerformAnalysis}


class PunchcardParser(object):  # Status: needs to be run but looks ok
    def __init__(self, root: str = "../punchcards") -> None:
        # NOTE: root should be changed to a directory defined in a ~/.cytograph_rc file
        self.root = root
        self._analyses_dict = {}  # type: Dict
        self.model = {}  # type: Dict
        self._load_model()
        self._load_defs()

    def _load_model(self) -> None:
        self.model = yaml.load(open(os.path.join(self.root, "Model.yaml")))

    def _load_defs(self) -> None:
        debug_msgs = defaultdict(list)  # type: dict
        for cur, dirs, files in os.walk(self.root):
            for file in files:
                if ((".yaml" in file) or (".yml" in file)) and ("Model.yaml" not in file):
                    temp_dict = yaml.load(open(os.path.join(self.root, file)))
                    name = temp_dict["abbreviation"]
                    model_copy = copy.deepcopy(self.model)

                    # Do an update of the model dictionary, so to keep the defaults
                    for k, v in self.model.items():
                        if type(v) == dict:
                            for kk, vv in v.items():
                                if type(vv) == dict:
                                    for kkk, vvv in vv.items():
                                        try:
                                            model_copy[k][kk][kkk] = temp_dict[k][kk][kkk]
                                        except KeyError:
                                            debug_msgs[name].append("Analysis %s `%s:%s:%s` was not found. The Default `%s` will be used" % (name, k, kk, kkk, model_copy[k][kk][kkk]))
                                else:
                                    try:
                                        model_copy[k][kk] = temp_dict[k][kk]
                                    except KeyError:
                                        debug_msgs[name].append("Analysis %s `%s:%s` was not found. The Default `%s` will be used" % (name, k, kk, model_copy[k][kk]))
                        else:
                            try:
                                model_copy[k] = temp_dict[k]
                            except KeyError:
                                debug_msgs[name].append("Analysis %s `%s` was not found. The Default `%s` will be used" % (name, k, model_copy[k]))
                    self._analyses_dict[name] = copy.deepcopy(model_copy)
                    self.debug_msgs = debug_msgs

    @property
    def all_analyses(self) -> List:
        return list(self._analyses_dict.values())

    @property
    def all_analyses_dict(self) -> Dict[str, Dict]:
        return dict(self._analyses_dict)

    def __getitem__(self, key: Any) -> Dict:
        for i in self.debug_msgs[key]:
            logging.debug(i)
        return self._analyses_dict[key]


def parse_punchcard_require(punchcard_obj: Dict) -> List[Tuple[luigi.Task]]:
    """Takes a dictionary parsed from the yaml file and returns the correspnding list of Tasks requirement

    The current version assumes that the input wil always be a WrapperTask
    """
    requirements: List[luigi.WrapperTask] = []
    for i in range(len(punchcard_obj["require"])):
        requirement_entry = punchcard_obj["require"][i]
        requirement_type = requirement_entry["type"]
        requirement_kwargs = requirement_entry["kwargs"]
        if requirement_type not in require_type_dict:
            raise NotImplementedError(f"Task type: {requirement_type} not allowed, you need to allow it adding it to require_type_dict")
        Task = require_type_dict[requirement_type]
        requirements += list(Task(**requirement_kwargs).requires())
    return requirements


def parse_punchcard_run(punchcard_obj: Dict) -> Iterator[luigi.Task]:
    """Yields luigi.Tasks after parsing out a dictionary describing the kind of tasks and their arguments
    """
    # the following safenames is implemented to make the gettattr statement secure
    safenames = set()  # type: set
    for k, v in dm.__dict__.items():
        if type(v) == luigi.task_register.Register:
            safenames |= {k}
    for task2run in punchcard_obj["run"]:
        task_type, task_kwargs = task2run["type"], task2run["kwargs"]
        if task_type not in safenames:
            raise NotImplementedError(f"Task type: {task_type} not allowed, becouse is not a valid luigi task")
        else:
            Task_class = getattr(dm, task_type)  # eval("cg.%s" % analysis_type)

            def Task(analysis: Any) -> luigi.Task:
                return Task_class(analysis, **task_kwargs)
            
            yield Task
