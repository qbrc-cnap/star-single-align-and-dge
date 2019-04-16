import json
import analysis.models
from analysis.models import AnalysisProject, ProjectConstraint

def check_constraints(project_pk, inputs_json_path):
    '''
    For this workflow, we can impose a constraint on the number of samples we will
    allow for processing.

    The length of the SingleEndRnaSeqAndDgeWorkflow.r1_files key will be used as the
    determinant of that number
    '''

    # load the inputs json:
    j = json.load(open(inputs_json_path))
    try:
        fastq_list = j['SingleEndRnaSeqAndDgeWorkflow.r1_files']
    except KeyError:
        # The chances of reaching this are very unlikely, but we are being extra careful here
        print('This should not happen-- the SingleEndRnaSeqAndDgeWorkflow.r1_files key should be present in your inputs JSON file')
        return False

    try:
        project = AnalysisProject.objects.get(pk=project_pk)
    except analysis.models.AnalysisProject.DoesNotExist:
        print('Query for analysis project with pk=%d was unsuccessful.  Something is quite wrong.' % project_pk)
        return False

    project_constraints = ProjectConstraint.objects.filter(project=project)
    constraint_num = len(project_constraints)
    # we should only have the single constraint for sample number:
    if constraint_num > 1:
        print('There should only be zero or one constraints on this type of workflow.  Please check what is going on')
        return False
    elif constraint_num == 0:
        return True
    else:
        # since we know that we imposed an AnalysisUnitConstraint, this is accessible as 'analysisunitconstraint' attribute
        # All classes deriving from ImplementedConstraint have the value attribute
        constraint_value = project_constraints[0].constraint.analysisunitconstraint.value

        # finally we can check if the constraints are satisfied:
        return len(fastq_list) <= constraint_value
