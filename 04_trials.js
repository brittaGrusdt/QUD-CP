// In this file you can specify the trial data for your experiment

// Test trials
var test_qud = `Based on the following dialogue between Ann and Bob, which of the two scenes is more likely to be the one that Bob is describing?` +
  `<br/><br/>Please select by clicking on the respective picture.<hr>`
var image_selection_trials = [];
TEST_DATA.forEach(function(dat) {
  let idx_quest = Object.values(QUESTS).indexOf(dat.question)
    image_selection_trials.push({
      id: dat.id,
      QUD: test_qud,
      question: Object.keys(QUESTS)[idx_quest],
      question_long: "<b>Ann</b>: " + dat.question,
      answer: "<b>Bob</b>: " + dat.answer,
      picture1: "stimuli/img/group/" + dat.id1 + ".png",
      picture2: "stimuli/img/group/" + dat.id2 + ".png",
      id1: dat.id1,
      id2: dat.id2,
      type: dat.type,
      causes_id1: dat.causes_id1,
      causes_id2: dat.causes_id2
    });
});
// adapt path to pictures depending on colour group in each trial
// add group and id separately
// here I changed code to save changes into slider_rating_trials 22.5. Malin
_.map(image_selection_trials, function (trial) {
  let group = _.sample(["group1", "group2"]);
  trial.picture1 = trial.picture1.replace("group", group);
  trial.picture2 = trial.picture2.replace("group", group);
  trial.group = group;
  trial.expected = TEST_EXPECT[trial.id]
});

EXAMPLE_TEST_TRIALS = {
  id: "example-test-trials",
  QUD: "EXAMPLE TRIAL<br/><br/>" + test_qud,
  picture1: "stimuli/img/example_test_trials1.png",
  picture2: "stimuli/img/example_test_trials2.png",
  id1: "example_test_trials1",
  id2: "example_test_trials2",
  answer: "<b>Bob</b>: " + ANSWERS.yellow + ANSWERS.normal_check,
  question_long: "<b>Ann</b>: " + QUESTS.neutral,
  question: "neutral",
  expected: "picture1",
  group: 'example',
  type: 'test-example',
  causes_id1: "",
  causes_id2: ""
}

// ----- TRAINING TRIALS (Buttons) for exp1 + exp2 ---- //
// the data of the training stimuli is always the same,
// the 4 buttons are always shown in same order
let TRAIN_TRIALS = [];
let train_keys = Object.keys(TRAIN_IDS);
train_keys.forEach(function (kind) {
  TRAIN_IDS[kind].forEach(function(id){
    let data = {
      QUD: "Use the buttons below to select all the blocks that you think will fall.<br/>(Multiple selections possible!)",
      id: id,
      question: '',
      expected: train_expectations[id],
      type: 'training'
    };
    TRAIN_TRIALS.push(data);
  });
});
