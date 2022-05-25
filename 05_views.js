// In this file you can instantiate your views
// We here first instantiate wrapping views, then the trial views

/** Wrapping views below

* Obligatory properties

    * trials: int - the number of trials this view will appear
    * name: string

*Optional properties
    * buttonText: string - the text on the button (default: 'next')
    * text: string - the text to be displayed in this view
    * title: string - the title of this view

    * More about the properties and functions of the wrapping views - https://magpie-ea.github.io/magpie-docs/01_designing_experiments/01_template_views/#wrapping-views

*/

// Every experiment should start with an intro view. Here you can welcome your participants and tell them what the experiment is about
const intro = magpieViews.view_generator("intro", {
  trials: 1,
  name: "intro",
  // If you use JavaScripts Template String `I am a Template String`, you can use HTML <></> and javascript ${} inside
  text: `Thank you for your participation in our study!
         Your anonymous data makes an important contribution to our understanding of human language use.
          <br />
          <br />
          Legal information:
          By answering the following questions, you are participating in a study
          being performed by scientists from the University of Osnabrueck.
          <br />
          <br />
          You must be at least 18 years old to participate.
          <br />
          <br />
          Your participation in this research is voluntary.
          You may decline further participation, at any time, without adverse consequences.
          <br />
          <br />
          Your anonymity is assured; the researchers who have requested your
          participation will not receive any personal information about you.
          `,
  buttonText: "begin the experiment"
});

const instructions_general = magpieViews.view_generator("instructions", {
  trials: 1,
  name: "instructions_general",
  title: "General Instructions",
  text: `This experiment consists of two phases, a <b>training</b> and a <b>testing</b> phase,
          with <b>8</b> respectively <b>18</b> trials.
         In total, you will need about <b>10 minutes</b> to finish.
         <br/>
         <br/>
         We will start with the training phase in which you will see arrangements of blocks.<br/>
         Your task will be to select all blocks that you think will fall (multiple selections possible!).
        <br/>
         After you have made your choice, click on RUN to see which blocks actually fall.
         <br/>
         By clicking on CONTINUE you'll proceed with the next training trial.
         <br/>
         <br/>
         <b>Please note</b>:
         <br/>
         We recommend to use <b>Full Screen Mode</b> throughout the experiment (usually switched on/off with F11).`,
  buttonText: "Start Training"
});

const instructions_test_text = magpieViews.view_generator("instructions", {
  trials: 1,
  name: "instructions_test_text",
  title: "Instructions Test Phase",
  text: `Great -  we'll now proceed with the test phase of the experiment.
    You will be shown two pictures of the same scenes of block arrangements as in the training phase.
    <br/>
    <br/>
    This time, you will read a dialogue between Ann and Bob talking about any such scene. Note that this scene is not necessarily one of the two pictures that you see.
    <br/>
    First, click on the upper button to see Ann's question, then click on the button below to see Bob's answer.
    <br />
    <br />
    After you read the dialogue, you will be able to select a picture by clicking on it.
    We ask you to select the picture that you think Bob is more likely to see given his answer to Ann's question.
    <br/>
    <br/>
    By clicking on CONTINUE you'll proceed with the next trial.
    <br />
    <br />
    There will be one example trial next before the actual test phase starts.
    `,
  buttonText: "go to example"
});

const instructions_test_img = magpieViews.view_generator("instructions", {
  trials: 1,
  name: "instructions_test_img",
  title: "Instructions Test Phase",
  text: `Great -  we'll now proceed with the test phase of the experiment.
    You will be shown two pictures of the same scenes of block arrangements as in the training phase and a dialog between Ann and Bob talking about such a scene.
    Note that this scene is not necessarily one of the two pictures that you see.
    <br/>
    <br/>
    Below, you see how test trials will look like. Follow the instructions given in 1. to 5. to complete a trial.
    <br />
    <br />
    <img src="./stimuli/img/instructions-test-anonym.jpeg" style="max-width:90%;">
    <br />
    <br />
    We includeded some simple attention check trials where you are explicitly told which picture to select.
    Note that if you do not follow these instructions, you will not get paid.
    <br />
    <br />
    There will be one example trial next before the actual test phase starts.
    `,
  buttonText: "go to example"
});


const post_test = magpieViews.view_generator("post_test", {
  trials: 1,
  name: "post_test",
  title: "Additional information",
  text: "Answering the following questions is optional, but your answers will help us analyze our results.",
  comments_question: 'Further comments'
});


// The 'thanks' view is crucial; never delete it; it submits the results!
const thanks = magpieViews.view_generator("thanks", {
  trials: 1,
  name: "thanks",
  title: "Thank you very much for taking part in this experiment!",
  prolificConfirmText: "Press the button"
});

const animation_view = animation_generator(
  {name: "animation_view",
   title: "title",
   CT: 0,
   trials: NB_TRAIN_TRIALS,
   data: ""}
 )

const test_image_selection = magpieViews.view_generator(
  "image_selection", {
    data: TEST_TRIALS,
    name: 'image_selection',
    trials: TEST_IDS.length
  }, {
    stimulus_container_generator: forced_choice_generator.stimulus_container_gen,
    answer_container_generator: forced_choice_generator.answer_container_gen,
    handle_response_function: forced_choice_generator.handle_response_function
  }
);

const example_image_selection = magpieViews.view_generator(
  "image_selection", {
    data: [EXAMPLE_TEST_TRIALS],
    name: 'image_selection',
    trials: 1
  }, {
    stimulus_container_generator: forced_choice_generator.stimulus_container_gen,
    answer_container_generator: forced_choice_generator.answer_container_gen,
    handle_response_function: forced_choice_generator.handle_response_function
  }
);

const divider_example_test_phase = magpieViews.view_generator("instructions", {
  trials: 1,
  name: "divider_example_test_phase",
  title: "Training completed!",
  text: `Please click on CONTINUE to start the actual experiment.`,
  buttonText: "Start Experiment"
});

// captcha
// speaker and listeneers names to be sampled from for the botcaptcha
var speaker = _.sample(["James", "John", "Robert", "Michael", "William", "David", "Richard", "Joseph", "Thomas", "Charles"]);
var listener = _.sample(["Mary", "Patricia", "Jennifer", "Linda", "Elizabeth", "Barbara", "Susan", "Jessica", "Sarah", "Margaret"]);

const botcaptcha = custom_botcaptcha({
  name: 'botcaptcha',
  trials: 1,
  story: speaker + ' says to ' + listener + ': "It\'s a beautiful day, isn\'t it?"',
  question: "Who is " + speaker + " talking to?",
  speaker: speaker,
  listener: listener
});
