// In this file you initialize and configure your experiment using magpieInit

$("document")
  .ready(function () {
    // prevent scrolling when space is pressed
    window.onkeydown = function (e) {
      if (e.keyCode === 32 && e.target === document.body) {
        e.preventDefault();
      }
    };

    // calls magpieInit
    // in debug mode this returns the magpie-object, which you can access in the console of your browser
    // e.g. >> window.magpie_monitor or window.magpie_monitor.findNextView()
    // in all other modes null will be returned
    window.magpie_monitor = magpieInit({
      // You have to specify all views you want to use in this experiment and the order of them
      //*alternating trials of Exp1 (prior elicitation) and Exp2 (production)* in views_test!//
      views_seq: [
      botcaptcha,
      intro,
      instructions_general,
      animation_view, // training phase
      instructions_test_img,
      example_image_selection, // test example
      divider_example_test_phase,
      test_image_selection, // test phase
      post_test,
      thanks
      ],
      // Here, you can specify all information for the deployment
      deploy: {
        experimentID: "67",
        serverAppURL: "https://mcmpact.ikw.uni-osnabrueck.de/magpie/api/submit_experiment/",
        // Possible deployment methods are:
        // "debug" and "directLink"
        // As well as "MTurk", "MTurkSandbox" and "Prolific"
        deployMethod: "Prolific",
        contact_email: "britta.grusdt@uni-osnabrueck.de",
        prolificURL: "https://app.prolific.co/submissions/complete?cc=38EEABF8"
      },
      // Here, you can specify how the progress bar should look like
      progress_bar: {
        // list the view-names of the views for which you want a progress bar
        in: [test_image_selection.name, animation_view.name],
      // Possible styles are "default", "separate" and "chunks"
        style: "default",
        width: 100
      }
    });
  });
