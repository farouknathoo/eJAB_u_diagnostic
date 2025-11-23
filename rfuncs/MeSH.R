# Step 2: Function to send one batch of up to 100 IDs
query_mesh_batch <- function(id_batch) {
  formatted_ids <- paste0("mesh:", id_batch, collapse = " ")
  query_string <- paste0(
    'PREFIX meshv: <http://id.nlm.nih.gov/mesh/vocab#> ',
    'PREFIX mesh: <http://id.nlm.nih.gov/mesh/> ',
    'SELECT ?meshID ?treeNumber WHERE { VALUES ?meshID { ',
    formatted_ids,
    ' } ?meshID meshv:treeNumber ?treeNumber . }'
  )

  query_url <- paste0(
    "https://id.nlm.nih.gov/mesh/sparql?query=",
    URLencode(query_string, reserved = TRUE),
    "&format=CSV&inference=false"
  )

  res <- GET(query_url)
  stop_for_status(res)

  csv_text <- content(res, as = "text", encoding = "UTF-8")
  read_csv(csv_text, show_col_types = FALSE)
}
