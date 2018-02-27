#include "moves.h"
//#include "read_probability_calculator.h"
#include <algorithm>
#include <cassert>
#include "graph.h"

int FindPathWithSameEnding(vector<Path>& paths, int pi) {
  for (size_t i = 0; i < paths.size(); i++) {
    auto &p = paths[i];
    auto rp = p.GetReverse();
    if ((int)i == pi) continue;
    if (p[0] == paths[pi].back()) {
      return i;
    }
    if (rp[0] == paths[pi].back()) {
      paths[i] = rp;
      return i;
    }
  }
  return -1;
}

void MergePaths(vector<Path>& paths, int pi, int same_end) {
  assert(paths[pi].back() == paths[same_end][0]);
  paths[pi].AppendPath(paths[same_end], 1);
  paths[pi].history_ += PATH_EXTEND_RANDOMLY;
  swap(paths[same_end], paths[paths.size()-1]);
  paths.pop_back();
}

bool ExtendPathsRandomly(const vector<Path>& paths, vector<Path>& out_paths,
                         const MoveConfig& config) {
  int pi = rand()%paths.size();
  out_paths = paths;
  random_shuffle(out_paths.begin(), out_paths.end());
  if (rand()%2 == 1) {
    out_paths[pi].Reverse();
  }
  if (!out_paths[pi].ExtendRandomly(config.big_node_threshold,
                                    config.rand_extend_step_threshold,
                                    config.rand_extend_distance_threshold)) {
    return false;
  }
  int same_end = FindPathWithSameEnding(out_paths, pi);
  if (same_end == -1) {
    return true;
  }

  MergePaths(out_paths, pi, same_end);
  return true;
}

bool BreakPaths(const vector<Path>& paths, vector<Path>& out_paths,
                const MoveConfig& config) {
  int pi = rand()%paths.size();
  out_paths = paths;
  if (out_paths[pi].size() < 2) {
    return false;
  }
  int break_pos = 1+rand()%(out_paths[pi].size()-1);
  Path p2 = out_paths[pi].CutAt(break_pos, config.big_node_threshold);
  out_paths[pi].history_ = out_paths[pi].history_ + PATH_BREAK;
  p2.history_ = out_paths[pi].history_;
  out_paths.push_back(p2);
  return true;
}

int chooseWeightedRandomly(const vector<int> &weights) {
  int total_weight = 0;
  for (int w: weights) total_weight += w;
  if (total_weight == 0) return -1;

  int pref_weight = 0;
  int res = -1;

  const int rw = rand() % total_weight;

  while (pref_weight <= rw) {
    res++;
    pref_weight += weights[res];
  }

  return res;
}

Path SampleRandomWalk(Node* start, const unordered_set<int>& allowed_nodes, const unordered_set<int>& desired_targets, const int MAX_LENGTH) {
  vector<Node*> nodes;
  nodes.push_back(start);
  int length = 0;

  vector<Node*> available;
  while (length < MAX_LENGTH && desired_targets.count(nodes.back()->id_) == 0) {
    Node* x = nodes.back();
    available.clear();
    for (auto nx: x->next_) {
      //cerr << nx->id_ << " ";
      if (allowed_nodes.count(nx->id_)) {
        available.push_back(nx);
      }
    }
    if (available.empty()) break;
    const int choice = rand() % available.size();
    //cerr << "choice: " << choice << endl;
    Node* chosen = available[choice];
    nodes.push_back(chosen);
    length += chosen->str_.size();
  }
  auto ret = Path(nodes);
  //cerr << "(UN?)POSSIBLE PATH: " << ret.ToDebugString() << endl;
  return ret;
}

bool JoinWithAdvicePaired(const vector<Path>& paths, vector<Path>& out_paths,
                    const MoveConfig& config, PairedReadProbabilityCalculator& pc) {
  cerr << "JoinWithAdvicePaired()" << endl;

  if (paths.empty()) return false;
  // choose path to extend (by joining another path) randomly
  // @TODO maybe take into consideration amount of mismatched paired reads
  const int start_path_id = rand() % (int)paths.size();
  const Path& start_path = paths[start_path_id];

  //cerr << "start_path_id: " << start_path_id << endl;
  //cerr << "start_path: " << start_path.ToDebugString() << endl;

  // exclude nondisjoint paths
  vector<int> disjoint_path_ids;
  const int disjoint_length = (int)(pc.mean_distance_ + pc.std_distance_ * 3);
  for (int i = 0; i < (int)paths.size(); i++) {
    if (start_path.isPartlyDisjoint(paths[i], disjoint_length)) {
      disjoint_path_ids.push_back(i);
    }
  };
  if (disjoint_path_ids.empty()) return false;

  //cerr << "Disjoint path ids:";
  //for (int i: disjoint_path_ids) cerr << " " << i;
  //cerr << endl;

  // align paired reads to the start path
  vector<SingleReadAlignment> start_left_alignments = pc.path_aligner_.GetPartAlignmentsForPath(start_path, 0);
  sort(start_left_alignments.begin(), start_left_alignments.end());
  vector<SingleReadAlignment> start_right_alignments = pc.path_aligner_.GetPartAlignmentsForPath(start_path, 1);
  sort(start_right_alignments.begin(), start_right_alignments.end());

  //cerr << "paired reads aligned to the start path" << endl;

  // align all paired reads to the filtered paths
  vector<vector<SingleReadAlignment>> left_alignments(paths.size()), right_alignments(paths.size());

  int counter = 0;
  for (int i: disjoint_path_ids) {
    left_alignments[i] = pc.path_aligner_.GetPartAlignmentsForPath(paths[i], 0);
    sort(left_alignments[i].begin(), left_alignments[i].end());
    right_alignments[i] = pc.path_aligner_.GetPartAlignmentsForPath(paths[i], 1);
    sort(right_alignments[i].begin(), right_alignments[i].end());

    // debug output
    counter++;
    if (counter % 10 == 0) cerr << "\r" << counter << " paths aligned";
  }
  cerr << endl; // debug output

  //cerr << "paired reads aligned to the filtered disjoint paths" << endl;

  // evaluate score for each path
  unordered_map<int,int> start_left_count, start_right_count;
  int max_read_id = 0;
  for (auto al: start_left_alignments){
    start_left_count[al.read_id] = 1 + start_left_count[al.read_id];
    max_read_id = max(max_read_id, al.read_id);
  }
  for (auto al: start_right_alignments) {
    start_right_count[al.read_id] = 1 + start_right_count[al.read_id];
    max_read_id = max(max_read_id, al.read_id);
  }

  //cerr << "alignment counts for start path evaluated." << endl;

  vector<int> path_score(paths.size(), 0);

  for (int path_id: disjoint_path_ids) {
    const auto &lefts = left_alignments[path_id];
    const auto &rights = right_alignments[path_id];
    if (lefts.empty() && rights.empty()) continue;

    auto it_left = lefts.begin();
    auto it_right = rights.begin();

    for (int read_id = 0; read_id <= max_read_id; read_id++) {
      // amount of left and right parts of read (id=read_id), aligned to the given (target) path
      int target_lefts = 0, target_rights = 0;

      while (it_left != lefts.end() && it_left->read_id < read_id) it_left++;
      while (it_left != lefts.end() && it_left->read_id == read_id) {
        target_lefts++;
        it_left++;
      }

      while (it_right != rights.end() && it_right->read_id < read_id) it_right++;
      while (it_right != rights.end() && it_right->read_id == read_id) {
        target_rights++;
        it_right++;
      }

      path_score[path_id] += start_left_count[read_id] * target_rights + start_right_count[read_id] * target_lefts;
    }
  }

  //cerr << "DISJOINT PATHS WITH SCORES:" << endl;
  //for (int path_id: disjoint_path_ids) {
  //  if (path_score[path_id] > 0)  cerr << path_score[path_id] << "\t!! " << paths[path_id].ToDebugString() << endl;
  //  //else cerr << "-\t   " << paths[path_id].ToDebugString() << endl;
  //}

  // choose randomly (based on score) the walk to join (with the orientation and complementarity)
  const int target_path_id = chooseWeightedRandomly(path_score);
  if (target_path_id == -1) return false;
  const Path& target_path = paths[target_path_id];
  //cerr << "Chosen target path: " << target_path_id << " :: " << target_path.ToDebugString() << endl;

  // sample (randomly) possible connections
  // choose the best

  // start path:  xb ----> xe. xbr := xb->rc_, xer := xe->rc_
  // target path: yb ----> ye. ybr := yb->rc_, yer := ye->rc_
  Node* xb = start_path.nodes_.front();
  Node* xe = start_path.nodes_.back();
  //Node* xbr = xb->rc_;
  Node* xer = xe->rc_;
  Node* yb = target_path.nodes_.front();
  Node* ye = target_path.nodes_.back();
  //Node* ybr = yb->rc_;
  Node* yer = ye->rc_;

  // 4 possibilities: xe ---> (yb | yer) ;  ye ---> ( xb | xer)
  Graph* G = xb->graph_;
  const vector<Node*> yb_drb = G->GetDrainageBasinForNode(yb);
  const vector<Node*> yer_drb = G->GetDrainageBasinForNode(yer);

  const vector<Node*> xb_drb = G->GetDrainageBasinForNode(xb);
  const vector<Node*> xer_drb = G->GetDrainageBasinForNode(xer);

  //cerr << "yb_drb\t:: "; for(auto x: yb_drb) cerr << x->id_ << " "; cerr << endl;
  //cerr << "ybr_drb\t:: "; for(auto x: ybr_drb) cerr << x->id_ << " "; cerr << endl;
  //cerr << "xb_drb\t:: "; for(auto x: xb_drb) cerr << x->id_ << " "; cerr << endl;
  //cerr << "xbr_drb\t:: "; for(auto x: xbr_drb) cerr << x->id_ << " "; cerr << endl;

  unordered_set<int> first_pool_ids, second_pool_ids;

  for (const auto &A: {yb_drb, yer_drb}) {
    if (find(A.begin(), A.end(), xe) != A.end()) {
      for (auto n: A) {
        first_pool_ids.insert(n->id_);
      }
    }
  }

  for (const auto &A: {xb_drb, xer_drb}) {
    if (find(A.begin(), A.end(), ye) != A.end()) {
      for (auto n: A) {
        second_pool_ids.insert(n->id_);
      }
    }
  }

  //cerr << "first pool ids: [size: " << first_pool_ids.size() << "] "; for(int x: first_pool_ids) cerr << x << " "; cerr << endl;
  //cerr << "second pool ids: [size: " << second_pool_ids.size() << "] "; for(int x: second_pool_ids) cerr << x << " "; cerr << endl;

  vector<Path> possible_connections;

  const int SAMPLE_NUM = 100;
  const int MAX_CONN_LENGTH = disjoint_length * 3; // bases, not nodes

  if (!first_pool_ids.empty()) {
    unordered_set<int> desired_targets({yb->id_, yer->id_});
    for (int num = 0; num < SAMPLE_NUM; num++) {
      Path p = SampleRandomWalk(xe, first_pool_ids, desired_targets, MAX_CONN_LENGTH);
      if (p.back() == yb || p.back() == yer) {
        bool is_new_path = true;
        for (auto ex_p: possible_connections) {
          if (ex_p.IsSameNoReverse(p)) {
            is_new_path = false;
            break;
          }
        }
        if (is_new_path) possible_connections.push_back(p);
      }
    }
  }
  if (!second_pool_ids.empty()) {
    unordered_set<int> desired_targets({xb->id_, xer->id_});
    for (int num = 0; num < SAMPLE_NUM; num++) {
      Path p = SampleRandomWalk(ye, second_pool_ids, desired_targets, MAX_CONN_LENGTH);
      if (p.back() == xb || p.back() == xer) {
        bool is_new_path = true;
        for (const auto ex_p: possible_connections) {
          if (ex_p.IsSameNoReverse(p)) {
            is_new_path = false;
            break;
          }
        }
        if (is_new_path) possible_connections.push_back(p);
      }
    }
  }

  //cerr << "possible connections: " << endl;
  //for (auto &p: possible_connections) {
  //  cerr << p.ToDebugString() << endl;
  //}

  if (possible_connections.empty()) return false;

  Path best_path;
  double best_score = -10000000000;

  for (auto &p: possible_connections) {
    vector<Node*> cut_nodes(p.nodes_.begin() + 1, p.nodes_.end() - 1);
    Path inter_path(cut_nodes);

    Path np;

    if (p[0] == xe) {
      if (p.nodes_.back() == yb) {
        np = Path(start_path.nodes_);
        np.AppendPath(inter_path);
        np.AppendPath(target_path);
      }
      if (p.nodes_.back() == yer) {
        np = Path(start_path.nodes_);
        np.AppendPath(inter_path);
        np.AppendPath(target_path.GetReverse());
      }
    }
    if (p[0] == ye) {
      if (p.nodes_.back() == xb) {
        np = Path(target_path.nodes_);
        np.AppendPath(inter_path);
        np.AppendPath(start_path);
      }
      if (p.nodes_.back() == xer) {
        np = Path(target_path.nodes_);
        np.AppendPath(inter_path);
        np.AppendPath(start_path.GetReverse());
      }
    }
    np.history_ = start_path.history_ + PATH_JOIN;

    // @TODO check for correct initialising of pp_change
    PairedProbabilityChange pp_change;
    pp_change.added_paths.push_back(np);
    pp_change.removed_paths.push_back(target_path);
    pp_change.removed_paths.push_back(start_path);

    pp_change.new_paths_length = pc.GetPathsLength(paths) - (int)start_path.ToString(true).size() - (int)target_path.ToString(true).size() + (int)np.ToString(true).size();
    //pp_change.new_paths = paths;

    pc.EvalProbabilityChange(pp_change, false);
    double score = pc.EvalTotalProbabilityFromChange(pp_change, false);
    cerr << "Path (score: "<< score << "): " << np.ToDebugString() << endl;
    if (score > best_score) {
      best_score = score;
      best_path = np;
    }
  }

  // add to the `out_paths`, return `true` value
  out_paths.clear();
  out_paths.reserve(paths.size() - 1);
  for (int i = 0; i < (int)paths.size(); i++) {
    if (i != start_path_id && i != target_path_id) {
      out_paths.push_back(paths[i]);
    }
  }
  out_paths.push_back(best_path);
  return true;
}

bool JoinWithAdvice(const vector<Path>& paths, vector<Path>& out_paths,
                    const MoveConfig& config, GlobalProbabilityCalculator& probability_calculator) {
  // choose dataset for advices. if no available, then return false
  vector<int> available_paired_read_sets;
  for (int i = 0; i < (int)probability_calculator.paired_read_calculators_.size(); i++) {
    auto &p = probability_calculator.paired_read_calculators_[i];
    if (p.first.use_as_advice_) {
      available_paired_read_sets.push_back(i);
    }
  }

  if (available_paired_read_sets.empty()) {
    return false;
  }

  const int which_dataset_to_choose = rand() % (int)available_paired_read_sets.size();
  return JoinWithAdvicePaired(paths, out_paths, config,
                              probability_calculator.paired_read_calculators_[available_paired_read_sets[which_dataset_to_choose]].first);

  // @TODO add HiC choice
}

void MakeMove(const vector<Path>& paths, vector<Path>& out_paths, const MoveConfig& config, GlobalProbabilityCalculator& probability_calculator,
              bool& accept_higher_prob) {
  while (true) {
    out_paths.clear();
    if (TryMove(paths, out_paths, config, probability_calculator, accept_higher_prob)) return;
  }
}

bool TryMove(const vector<Path>& paths, vector<Path>& out_paths, const MoveConfig& config, GlobalProbabilityCalculator& probability_calculator,
             bool& accept_higher_prob) {
  // @TODO add probs of moves into config

  int move = rand()%5;
  if (move < 2) {
    accept_higher_prob = false;
    return ExtendPathsRandomly(paths, out_paths, config);
  }
  if (2 <= move && move < 4) {
    accept_higher_prob = true;
    return BreakPaths(paths, out_paths, config);
  }
  // @TODO Joining with advice move (high priority)
  if (move == 4) {
    accept_higher_prob = false; // @TODO check with Usama if correct
    return JoinWithAdvice(paths, out_paths, config, probability_calculator);
  }

  // @TODO add Local improvement move (low priority)
  // @TODO add Repeat optimization move (low priority)
  // @TODO Repeat interchange move (low priority)
  return false;
}
