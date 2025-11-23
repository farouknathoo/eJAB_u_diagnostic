import json
import csv
import hashlib

def generate_unique_id(nctId, analysis, outcome):
    """
    Generate a unique hash-based ID for each analysis based on study ID, entire outcome data,
    and entire analysis data.
    """
    # Convert the entire outcome data to a JSON string for consistent hashing
    outcome_json_str = json.dumps(outcome, sort_keys=True, ensure_ascii=False)
    analysis_json_str = json.dumps(analysis, sort_keys=True, ensure_ascii=False)
    unique_string = f"{nctId}_{outcome_json_str}_{analysis_json_str}"
    return hashlib.md5(unique_string.encode()).hexdigest()

def flatten_record(data):
    """
    Flatten a single JSON record to create a new row for each combination of:
      - outcome measure (only if it has analyses)
      - each analysis (pValue and statisticalMethod)
      - each groupId/count pair from denominators

    Also includes MESH terms for conditions and interventions, along with studyType,
    one-hot encoded phases, and additional designInfo fields: interventionModel,
    observationalModel, timePerspective, and allocation.
    """
    records = []
    
    # Extract basic study information
    protocol = data.get("protocolSection", {})
    identification = protocol.get("identificationModule", {})
    design = protocol.get("designModule", {})
    nctId = identification.get("nctId")
    briefTitle = identification.get("briefTitle")
    
    # Extract studyType and phases from designModule
    study_type = design.get("studyType", "")
    phases = design.get("phases", [])  # Expecting a list of phase strings
    # Define the set of phases to one-hot encode
    possible_phases = ["EARLY_PHASE1","PHASE1", "PHASE2", "PHASE3", "PHASE4"]
    
    # Extract additional designInfo fields (if available) from designModule
    design_info = design.get("designInfo", {})
    interventionModel = design_info.get("interventionModel", "")
    observationalModel = design_info.get("observationalModel", "")
    timePerspective = design_info.get("timePerspective", "")
    allocation = design_info.get("allocation", "")
    
    # Extract MESH terms for conditions and interventions
    derived = data.get("derivedSection", {})
    condition_meshes = derived.get("conditionBrowseModule", {}).get("meshes", [])
    condition_terms = "; ".join([mesh.get("term", "") for mesh in condition_meshes if mesh.get("term")])
    condition_ids   = "; ".join([mesh.get("id", "") for mesh in condition_meshes if mesh.get("id")])

    intervention_meshes = derived.get("interventionBrowseModule", {}).get("meshes", [])
    intervention_terms = "; ".join([mesh.get("term", "") for mesh in intervention_meshes if mesh.get("term")])
    intervention_ids   = "; ".join([mesh.get("id", "") for mesh in intervention_meshes if mesh.get("id")])

    # Determine if the study has at least one reference of type "RESULTS" or "DERIVED"
    references = protocol.get("referencesModule", {}).get("references", [])
    hasResultsOrDerived = any(ref.get("type") in ["RESULTS", "DERIVED"] for ref in references)
    
    # Iterate over outcome measures
    outcome_measures = data.get("resultsSection", {}).get("outcomeMeasuresModule", {}).get("outcomeMeasures", [])
    
    for outcome in outcome_measures:
        # Skip outcome measure if no analyses are present
        analyses = outcome.get("analyses")
        if not analyses:
            continue
        
        outcome_type = outcome.get("type")
        outcome_title = outcome.get("title")
        outcome_description = outcome.get("description")
        
        for analysis in analyses:
            p_value = analysis.get("pValue", "")
            stat_method = analysis.get("statisticalMethod", "")
            group_description = analysis.get("groupDescription", "")
            group_ids = analysis.get("groupIds", [])
            # Generate a unique ID for this analysis
            analysis_id = generate_unique_id(nctId, analysis, outcome)
            
            denoms = outcome.get("denoms", [])
            for denom in denoms:
                units = denom.get("units")
                if units != "Participants":
                    continue
                counts = denom.get("counts", [])
                for count in counts:
                    groupId = count.get("groupId")
                    # Only include this groupId if it's in the analysis's groupIds list
                    if group_ids and groupId not in group_ids:
                        continue
                        
                    value = count.get("value")
                    
                    # Create the record dictionary with the extra designInfo fields
                    record = {
                        "analysisId": analysis_id,
                        "nctId": nctId,
                        "briefTitle": briefTitle,
                        "outcomeType": outcome_type,
                        "outcomeTitle": outcome_title,
                        "units": units,
                        "groupId": groupId,
                        "value": value,
                        "pValue": p_value,
                        "statisticalMethod": stat_method,
                        "condition_mesh": condition_terms,
                        "condition_ids": condition_ids,
                        "intervention_mesh": intervention_terms,
                        "intervention_ids": intervention_ids,
                        "hasResultsOrDerived": hasResultsOrDerived,
                        "studyType": study_type,
                        "interventionModel": interventionModel,
                        "observationalModel": observationalModel,
                        "timePerspective": timePerspective,
                        "allocation": allocation,
                        "groupDescription": group_description
                    }
                    # One-hot encode the phases: add a new key for each possible phase
                    for phase in possible_phases:
                        record[phase] = 1 if phase in phases else 0
                    records.append(record)
    
    return records

# Read JSON from file
with open("ctg-studies_082025.json", "r", encoding="utf-8") as infile:
    json_data = json.load(infile)

# Process and flatten the JSON records
flattened_records = []
if isinstance(json_data, list):
    for record in json_data:
        flattened_records.extend(flatten_record(record))
else:
    flattened_records = flatten_record(json_data)

# Write the flattened records to a CSV file
with open("study3.csv", "w", newline="", encoding="utf-8") as csvfile:
    fieldnames = [
        "analysisId", "nctId", "briefTitle", "outcomeType", "outcomeTitle", "groupDescription",
        "units", "groupId", "value", "pValue", "statisticalMethod",
        "condition_mesh", "condition_ids","intervention_mesh", "intervention_ids", "hasResultsOrDerived",
        "studyType", "interventionModel", "observationalModel",
        "timePerspective", "allocation", "EARLY_PHASE1",
        "PHASE1", "PHASE2", "PHASE3", "PHASE4"
    ]
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    for rec in flattened_records:
        writer.writerow(rec)

print("CSV conversion complete. See study3.csv for results.")
