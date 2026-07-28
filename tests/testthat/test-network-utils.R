## Create a fake MetaPathNet-style network (no KEGG access needed)
fake_net <- matrix(
  c("cpd:C00078", "K00453", "k_compound:irreversible",
    "K00453", "cpd:C00328", "k_compound:irreversible",
    "cpd:C00328", "K00486", "k_compound:reversible",
    "K00486", "cpd:C00328", "k_compound:reversible",
    "K00486", "cpd:C00780", "k_compound:irreversible",
    "dr:D00001", "K00453", "k_compound:irreversible"),
  ncol = 3, byrow = TRUE
)
colnames(fake_net) <- c("source", "target", "interaction_type")

## -----------------------------------------------------------
## findMappedNodes
## -----------------------------------------------------------
test_that("findMappedNodes finds mapped and unmapped nodes", {
  res <- findMappedNodes(
    nodes = c("cpd:C00078", "K00453", "FAKE_NODE"),
    network_table = fake_net
  )

  expect_true("cpd:C00078" %in% res$mapped_nodes)
  expect_true("K00453" %in% res$mapped_nodes)
  expect_true("FAKE_NODE" %in% res$unmapped_nodes)
})

test_that("findMappedNodes normalises Cxxxxx to cpd:Cxxxxx", {
  res <- findMappedNodes(
    nodes = c("C00078"),
    network_table = fake_net
  )

  expect_true("cpd:C00078" %in% res$mapped_nodes)
})

## -----------------------------------------------------------
## combineNetworks
## -----------------------------------------------------------
test_that("combineNetworks merges two networks", {
  net_a <- fake_net[1:3, , drop = FALSE]
  net_b <- fake_net[4:6, , drop = FALSE]

  merged <- combineNetworks(net_a, net_b)

  expect_true(is.matrix(merged))
  expect_equal(ncol(merged), 3)
  expect_true(nrow(merged) >= max(nrow(net_a), nrow(net_b)))
})

test_that("combineNetworks removes duplicate edges", {
  merged <- combineNetworks(fake_net, fake_net)

  expect_equal(nrow(merged), nrow(unique(fake_net)))
})

test_that("combineNetworks rejects single input", {
  expect_error(combineNetworks(fake_net))
})

## -----------------------------------------------------------
## removeNetworkNodes
## -----------------------------------------------------------
test_that("removeNetworkNodes removes a node and its edges", {
  result <- removeNetworkNodes(
    nodes_to_remove = "K00453",
    network_table   = fake_net
  )

  all_nodes <- unique(as.vector(result[, 1:2]))
  expect_false("K00453" %in% all_nodes)
  expect_true(nrow(result) < nrow(fake_net))
})

test_that("removeNetworkNodes normalises compound IDs", {
  result <- removeNetworkNodes(
    nodes_to_remove = "C00078",
    network_table   = fake_net
  )

  all_nodes <- unique(as.vector(result[, 1:2]))
  expect_false("cpd:C00078" %in% all_nodes)
})

## -----------------------------------------------------------
## removeDrugs
## -----------------------------------------------------------
test_that("removeDrugs removes dr: nodes", {
  result <- removeDrugs(fake_net)

  all_nodes <- unique(as.vector(result[, 1:2]))
  has_drug <- any(grepl("^dr:", all_nodes))
  expect_false(has_drug)
  expect_true(nrow(result) < nrow(fake_net))
})

test_that("removeDrugs returns unchanged network when no drugs present", {
  no_drug_net <- fake_net[!grepl("dr:", fake_net[, 1]) & !grepl("dr:", fake_net[, 2]), , drop = FALSE]
  result <- removeDrugs(no_drug_net)

  expect_equal(nrow(result), nrow(no_drug_net))
})

## -----------------------------------------------------------
## replaceNode
## -----------------------------------------------------------
test_that("replaceNode replaces a node ID", {
  result <- replaceNode(
    nodes_to_replace = "K00453",
    replacement_node = "K99999",
    network_table    = fake_net
  )

  all_nodes <- unique(as.vector(result[, 1:2]))
  expect_false("K00453" %in% all_nodes)
  expect_true("K99999" %in% all_nodes)
})

test_that("replaceNode normalises compound IDs", {
  result <- replaceNode(
    nodes_to_replace = "C00078",
    replacement_node = "C99999",
    network_table    = fake_net
  )

  all_nodes <- unique(as.vector(result[, 1:2]))
  expect_false("cpd:C00078" %in% all_nodes)
  expect_true("cpd:C99999" %in% all_nodes)
})
