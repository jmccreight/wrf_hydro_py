import ast
from boltons.iterutils import remap, get_path
import copy
#import datetime
#import multiprocessing
import pathlib
from typing import Union
#import pickle
#import shlex
#import subprocess
#import time
#import uuid

from .ensemble_tools import DeepDiffEq, dictify, get_sub_objs

from .job import Job
from .schedulers import Scheduler
from .simulation import Simulation

#from .job_tools import solve_model_start_end_times

# Classes for constructing and running a wrf_hydro simulation
class EnsembleSimulation(object):
    """ TODO
    """

    def __init__(self):
        """ TODO """

        self.members = []
        """list: a list of simulations which are the members of the ensemble."""

        self.__member_diffs = {}
        """dict: a dictionary containing the differences across all the members attributes."""

        self.jobs = []
        """list: a list containing Job objects"""

        self.scheduler = None
        """Scheduler: A scheduler object to use for each Job in self.jobs"""

    def __len__(self):
        return(len(self.members))

    # The "canonical" name for len
    @property
    def N(self):
        return(self.__len__())

    # Metadata to store with the "member" simulations, conceptually this
    # data belongs to the members:
    # 1) member number
    # 2) description
    # 3) member_dir
    # 4) forcing_source_dir

    def add(
        self,
        obj: Union[list, Scheduler, Job]
    ):
        """Add an approparite object to an EnsembleSimulation, such as a Simulation, Job, or 
        Scheduler"""
        if isinstance(obj, list) or isinstance(obj, Simulation):
            self._addsimulation(obj)
        elif issubclass(type(obj), Scheduler):
            self._addscheduler(obj)
        elif isinstance(obj, Job):
            self._addjob(obj)
        else:
            raise TypeError('obj is not of a type expected for a EnsembleSimulation')

    def _addscheduler(self, scheduler: Scheduler):
        """Private method to add a Scheduler to a Simulation
        Args:
            scheduler: The Scheduler to add
        """
        self.scheduler = copy.deepcopy(scheduler)

    def _addjob(self, job: Job):
        """Private method to add a job to a Simulation
        Args:
            job: The job to add
        """
        job = copy.deepcopy(job)
        # Postpone the following until compose and do it on the
        # individual members.
        # job._add_hydro_namelist(self.base_hydro_namelist)
        # job._add_hrldas_namelist(self.base_hrldas_namelist)
        self.jobs.append(job)


    def _addsimulation(
        self,
        sims: Union[list, Simulation]
    ):
        """Private method to add a Model to a Simulation
        Args:
            model: The Model to add
        """

        if(type(sims) is Simulation):
            sims = [copy.deepcopy(sims)]

        for nn in sims:

            if type(nn) is not Simulation:
                raise valueError("A non-simulation object can not be "
                                 "added to the ensemble members")

            self.members.append(copy.deepcopy(nn))

            # If copying an existing ensemble member, nuke the following metadatas
            nuke_metadata = ['number', 'jobs', 'scheduler']
            # number is the detector for all ensemble metadata.
            for metadatum in nuke_metadata:
                if hasattr(nn, metadatum):
                    delattr(self.members[len(self.members)-1], metadatum)

        # Put refs to these properties in the ensemble objects
        for mm in range(len(self.members)):
            if not hasattr(self.members[mm], 'number'):
                self.members[mm].number = "%03d" % (mm,)
                self.members[mm].run_dir = 'member_' + self.members[mm].number

    # A quick way to setup a basic ensemble from a single sim.
    def replicate_member(
        self,
        N: int,
        copy_members: bool=True
    ):
        if self.N > 1:
            raise ValueError('The ensemble must only have one member to replicate.')
        else:
            for nn in range(1,N):
                self.add(self.members[0])

    # -------------------------------------------------------
    # The member_diffs attribute has getter (@property) and setter methods.
    # The get method summarizes all differences across all the attributes of the
    #   members list attribute and (should) only report member attributes when there
    #   is at least one difference between members.
    # The setter method is meant as a convenient way to specify the differences in
    #   member attributes across the ensemble.

    @property
    def member_diffs(self):

        if len(self) == 1:
            print('Ensemble is of lenght 1, no differences.')
            return {}

        mem_0_ref_dict = dictify(self.members[0])

        all_diff_keys = set({})
        for ii in range(1, len(self)):
            mem_ii_ref_dict = dictify(self.members[ii])
            diff = DeepDiffEq(mem_0_ref_dict, mem_ii_ref_dict, eq_types={pathlib.PosixPath})

            unexpected_diffs = set(diff.keys()) - set(['values_changed'])
            if len(unexpected_diffs):
                unexpected_diffs1 = {uu: diff[uu] for uu in list(unexpected_diffs)}
                raise ValueError(
                    'Unexpected attribute differences between ensemble members:',
                    unexpected_diffs1
                )

            diff_keys = list(diff['values_changed'].keys())
            all_diff_keys = all_diff_keys | set([ss.replace('root', '') for ss in diff_keys])

        diff_tuples = [ss.replace('][', ',') for ss in list(all_diff_keys)]
        diff_tuples = [ss.replace('[', '(') for ss in list(diff_tuples)]
        diff_tuples = [ss.replace(']', ')') for ss in list(diff_tuples)]
        diff_tuples = [ast.literal_eval(ss) for ss in list(diff_tuples)]

        self.__member_diffs = {}
        for dd in diff_tuples:
            self.__member_diffs[dd] = [get_path(dictify(mm), dd) for mm in self.members]

        return(self.__member_diffs)

    def set_member_diffs(
        self,
        att_tuple: tuple,
        values: list
    ):

        if type(values) is not list:
            values = [values]

        if len(values) == 1:
            the_value = values[0]
            values = [the_value for ii in range(len(self))]

        if len(values) != len(self):
            raise ValueError("The number of values supplied does not equal the number of members.")

        def update_obj_dict(obj, att_tuple):

            def visit(path, key, value):
                superpath = path + (key,)

                if superpath != att_tuple[0:len(superpath)]:
                    return True
                if len(superpath) == len(att_tuple):
                    return key, new_value
                return True

            the_remap = remap(obj.__dict__, visit)
            obj.__dict__.update(the_remap)
            for ss in get_sub_objs(obj.__dict__):
                att_tuple_0 = att_tuple
                att_tuple = att_tuple[1:]
                if len(att_tuple) > 0:
                    update_obj_dict(obj.__dict__[ss], att_tuple)
                att_tuple = att_tuple_0

        for ii in range(len(self)):
            new_value = values[ii]
            update_obj_dict(self.members[ii], att_tuple)


# class WrfHydroEnsembleRun(object):
#     def __init__(
#         self,
#         ens_setup: WrfHydroEnsembleSetup,
#         run_dir: str,
#         rm_existing_run_dir = False,
#         mode: str='r',
#         jobs: list=None
#     ):

#         self.ens_setup = copy.deepcopy(ens_setup)
#         """WrfHydroSetup: The WrfHydroSetup object used for the run"""

#         # TODO(JLM): check all the setup members have to have rundirs with same path as run_dir
#         self.run_dir = pathlib.PosixPath(run_dir)
#         """pathlib.PosixPath: The location of where the jobs will be executed."""

#         self.jobs_completed = []
#         """Job: A list of previously executed jobs for this run."""
#         self.jobs_pending = []
#         """Job: A list of jobs *scheduled* to be executed for this run
#             with prior job dependence."""
#         self.job_active = None
#         """Job: The job currently executing."""
#         self.object_id = None
#         """str: A unique id to join object to run directory."""
#         self.members = []
#         """List: ensemble of Run Objects."""

#         # #################################

#         # Make run_dir directory if it does not exist.
#         if self.run_dir.is_dir() and not rm_existing_run_dir:
#             raise ValueError("Run directory already exists and rm_existing_run_dir is False.")

#         if self.run_dir.exists():
#             shutil.rmtree(str(self.run_dir))
#             self.run_dir.mkdir(parents=True)

#         # This sets up the runs. Writes WrfHydroRun.pkl objects to each dir.
#         for mm in self.ens_setup.members:
#             self.members.append(WrfHydroRun(
#                 mm,
#                 run_dir = mm.run_dir,
#                 deepcopy_setup=False
#             )
#         )

#         if jobs:
#             self.add_jobs(jobs)
#         else:
#             self.collect_output()
#             self.pickle()


#     def add_jobs(
#         self,
#         jobs: list
#     ):
#         """Add an Ensemble Run Job (array)."""

#         # Dont tamper with the passed object, let it remain a template in the calling level.
#         jobs = copy.deepcopy(jobs)

#         if type(jobs) is not list:
#             jobs = [jobs]

#         for jj in jobs:

#             # Attempt to add the job
#             if jj.scheduler:

#                 # A scheduled job can be appended to the jobs.pending list if
#                 # 1) there are no active or pending jobs
#                 # 2) if it is (made) dependent on the last active or pending job.

#                 # Get the job id of the last active or pending job.
#                 last_job_id = None
#                 if self.job_active:
#                     last_job_id = self.job_active.sched_job_id
#                 if len(self.jobs_pending):
#                     last_job_id = self.jobs_pending[-1].scheduler.sched_job_id

#                 # Check the dependency on a previous job
#                 if last_job_id is not None:
#                     if jj.scheduler.afterok is not None and jj.scheduler.afterok != last_job_id:
#                         raise ValueError("The job's dependency/afterok conflicts with reality.")
#                     jj.scheduler.afterok = last_job_id
#                 #else: 
#                 #    if jj.scheduler.afterok is not None:
#                 #        raise ValueError("The job's dependency/afterok conflicts with reality.")

#             # Set submission-time job variables here.
#             jj.user = get_user()
#             job_submission_time = datetime.datetime.now()
#             jj.job_submission_time = str(job_submission_time)
#             jj.job_date_id = 'job_' + str(len(self.jobs_completed) +
#                                           bool(self.job_active) +
#                                           len(self.jobs_pending))
#             # alternative" '{date:%Y-%m-%d-%H-%M-%S-%f}'.format(date=job_submission_time)
#             if jj.scheduler:
#                 jj.scheduler.array_size = len(self.members)

#             for mm in self.members:
#                 mm.add_jobs(jj)

#             self.jobs_pending.append(jj)

#         self.collect_output()
#         self.pickle()


#     def run_jobs(
#         self,
#         hold: bool=False,
#         n_mem_simultaneous: int=1
#     ):

#         hold_all = hold
#         del hold
        
#         # make sure all jobs are either scheduled or interactive?
                
#         if self.job_active is not None:
#             raise ValueError("There is an active ensemble run.")

#         if self.jobs_pending[0].scheduler:

#             run_dir = self.run_dir

#             # submit the jobs_pending.
#             job_afterok = None
#             hold = True

#             # For each the job arrays,
#             for ii, _ in enumerate(self.jobs_pending):

#                 #  For all the members,
#                 for mm in self.members:

#                     # Set the dependence into all the members jobs,
#                     jj = mm.jobs_pending[ii]
#                     jj.scheduler.afterok = job_afterok

#                     # Write everything except the submission script,
#                     # (Job has is_job_array == TRUE)
#                     jj.schedule(mm.run_dir, hold=hold)


#                 # Submit the array job for all the members, using the last member [-1] to do so.
#                 job_afterok = jj.schedule(self.run_dir, hold=hold, submit_array=True)
#                 # Keep that info in the object.
#                 self.jobs_pending[ii].sched_job_id = job_afterok
#                 # This is the "sweeper" job for job arrays.
#                 #job_afterok = self.jobs_pending[ii].collect_job_array(str_to_collect_ensemble)
#                 # Only hold the first job array
#                 hold = False

#             self.job_active = self.jobs_pending.pop(0)
#             self.pickle()

#             if not hold_all:
#                 self.members[-1].jobs_pending[0].release() # This prints the jobid
#             else:
#                 print(self.members[-1].jobs_pending[0].scheduler.sched_job_id)

#             self.destruct()
#             return run_dir

#         else:

#             # Make an attribute of this.
#             self.n_mem_simultaneous = n_mem_simultaneous
            
#             for jj in range(0, len(self.jobs_pending)):

#                 self.job_active = self.jobs_pending.pop(0)

#                 # This is a parallel a for loop over all members
#                 print("n_mem_simultaneous: ", n_mem_simultaneous, flush=True)
                
#                 pool = multiprocessing.Pool(n_mem_simultaneous)
#                 _ = pool.map(parallel_run_jobs, (mm for mm in self.members))
                    
#                 self.jobs_completed.append(self.job_active)
#                 self.job_active = None
                
#                 _ = pool.map(parallel_collect_jobs, (mm for mm in self.members))
#                 #self.collect_output()
#                 self.pickle()


#     str_to_collect_ensemble = (
#         "import pickle \n"
#         "import sys \n"
#         "import wrfhydropy \n"
#         "ens_run = pickle.load(open('WrfHydroEnsembleRun.pkl', 'rb')) \n",
#         "ens_run.collect_ensemble_runs(35) \n"
#         "sys.exit()"
#     )[0]


#     def collect_ensemble_runs(
#         self,
#         n_mem_simultaneous: int=1
#     ):
#         """Collect a completed job array. """

#         def n_jobs_not_complete(run_dir):
#             the_cmd = '/bin/bash -c "ls member_*/.job_not_complete 2> /dev/null | wc -l"'
#             stdout = subprocess.run(shlex.split(the_cmd), cwd=run_dir, stdout=subprocess.PIPE).stdout
#             ret = int(stdout.splitlines()[0].decode("utf-8"))
#             return ret

#         while n_jobs_not_complete(self.run_dir) != 0:
#             _ = time.sleep(.1)

#         if self.job_active:
#             self.jobs_completed.append(self.job_active)
#             self.job_active = None

#         #print('collect:',n_mem_simultaneous)
#         # Fairly minor difference between the speed with 80 members... 
#         pool = multiprocessing.Pool(n_mem_simultaneous)
#         self.members = pool.map(parallel_collectpickled_runs, (mm for mm in self.members))
#         #for ii, _ in enumerate(self.members):
#         #    self.members[ii] = self.members[ii].unpickle()
        
#         self.collect_output()
#         self.pickle()


#     def collect_output(self):
#         for mm in self.members:
#             mm.collect_output()


#     def pickle(self):
#         """Pickle the Run object into its run directory. Collect output first."""

#         # create a UID for the run and save in file
#         self.object_id = str(uuid.uuid4())
#         with open(self.run_dir.joinpath('.uid'), 'w') as f:
#             f.write(self.object_id)

#         # Save object to run directory
#         # Save the object out to the compile directory
#         with open(self.run_dir.joinpath('WrfHydroEnsembleRun.pkl'), 'wb') as f:
#             pickle.dump(self, f, 2)


#     def unpickle(self):
#         """ Load run object from run directory after a scheduler job. """
#         with open(self.run_dir.joinpath('WrfHydroEnsembleRun.pkl'), 'rb') as f:
#             return(pickle.load(f))


#     def destruct(self):
#         """ Pickle first. This gets rid of everything but the methods."""
#         self.pickle()
#         print("Jobs have been submitted to  the scheduler: This run object will now self destruct.")
#         self.__dict__ = {}
